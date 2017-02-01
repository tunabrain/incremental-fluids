/*
Copyright (c) 2013 Benedikt Bitterli

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.
*/

#define _USE_MATH_DEFINES

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <stack>

#include "../lodepng/lodepng.h"

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> int nsgn(T val) {
    return (val < T(0) ? -1 : 1);
}

double length(double x, double y) {
    return sqrt(x*x + y*y);
}

double cubicPulse(double x) {
    x = min(fabs(x), 1.0);
    return 1.0 - x*x*(3.0 - 2.0*x);
}

void rotate(double &x, double &y, double phi) {
    double tmpX = x, tmpY = y;
    x =  cos(phi)*tmpX + sin(phi)*tmpY;
    y = -sin(phi)*tmpX + cos(phi)*tmpY;
}

/* For three corners in a 1x1 square, with `in' being adjacent to `out1' and
 * `out2' and all three parameters being distances to a surface, `in' being
 * inside the surface and `out1' and `out2' outside, returns the area of the
 * square occupied by the surface.
 */
double triangleOccupancy(double out1, double in, double out2) {
	return 0.5*in*in/((out1 - in)*(out2 - in));
}

/* For four corners in a 1x1 square, with all parameters being distances to a
 * surface and `in1' and `in2 inside the surface, returns the are of the square
 * occupied by the surface.
 */
double trapezoidOccupancy(double out1, double out2, double in1, double in2) {
	return 0.5*(-in1/(out1 - in1) - in2/(out2 - in2));
}

/* Given the distance of four corners in a 1x1 square to a surface, returns the
 * area of the part of the square occupied by the surface computed analytically.
 *
 * The basic workings of this algorithm are quite similar to marching squares
 * (2D marching cubes). First, a mask is computed based on which corners are
 * inside and which are outside.
 * Based on this mask, the function differentiates between one of four cases:
 * a) Only one corner is inside
 *   => Compute using triangle area
 * b) Only one corner is outside
 *   => Invert distance field, compute 1 - triangle area
 * c) Two adjacent corners are inside
 *   => Compute using trapezoid area
 * d) Two opposing corners are inside
 *   => Compute as sum of area of two opposed triangles
 *
 * The two remaining cases, all corners outside/inside, can be computed trivially
 */
double occupancy(double d11, double d12, double d21, double d22) {
	double ds[] = {d11, d12, d22, d21};

    /* Compute mask */
	uint8_t b = 0;
	for (int i = 3; i >= 0; i--)
		b = (b << 1) | (ds[i] < 0.0 ? 1 : 0);
    
	switch (b) {
    /* All outside */
	case 0x0: return 0.0;
    /* One inside */
	case 0x1: return triangleOccupancy(d21, d11, d12);
	case 0x2: return triangleOccupancy(d11, d12, d22);
	case 0x4: return triangleOccupancy(d12, d22, d21);
	case 0x8: return triangleOccupancy(d22, d21, d11);
    /* One outside */
	case 0xE: return 1.0 - triangleOccupancy(-d21, -d11, -d12);
	case 0xD: return 1.0 - triangleOccupancy(-d11, -d12, -d22);
	case 0xB: return 1.0 - triangleOccupancy(-d12, -d22, -d21);
	case 0x7: return 1.0 - triangleOccupancy(-d22, -d21, -d11);
    /* Two adjacent inside */
	case 0x3: return trapezoidOccupancy(d21, d22, d11, d12);
	case 0x6: return trapezoidOccupancy(d11, d21, d12, d22);
	case 0x9: return trapezoidOccupancy(d12, d22, d11, d21);
	case 0xC: return trapezoidOccupancy(d11, d12, d21, d22);
    /* Two opposed inside */
	case 0x5: return triangleOccupancy(d11, d12, d22) +
		             triangleOccupancy(d22, d21, d11);
	case 0xA: return triangleOccupancy(d21, d11, d12) +
		             triangleOccupancy(d12, d22, d21);
    /* All inside */
	case 0xF: return 1.0;
	}
    
    return 0.0;
}

enum CellType {
    CELL_FLUID,
    CELL_SOLID
};

class SolidBody {
protected:
    double _posX;
    double _posY;
    double _scaleX;
    double _scaleY;
    double _theta;
    
    double _velX;
    double _velY;
    double _velTheta;
    
    void globalToLocal(double &x, double &y) const {
        x -= _posX;
        y -= _posY;
        rotate(x, y, -_theta);
        x /= _scaleX;
        y /= _scaleY;
    }
    
    void localToGlobal(double &x, double &y) const {
        x *= _scaleX;
        y *= _scaleY;
        rotate(x, y, _theta);
        x += _posX;
        y += _posY;
    }
    
    SolidBody(double posX, double posY, double scaleX, double scaleY,
        double theta, double velX, double velY, double velTheta) :
            _posX(posX), _posY(posY), _scaleX(scaleX), _scaleY(scaleY),
            _theta(theta), _velX(velX), _velY(velY), _velTheta(velTheta) {}
                
    virtual ~SolidBody() {};
    
public:
    virtual double distance(double x, double y) const = 0;
    virtual void closestSurfacePoint(double &x, double &y) const = 0;
    virtual void distanceNormal(double &nx, double &ny, double x, double y) const = 0;
    
    double velocityX(double x, double y) const {
        return (_posY - y)*_velTheta + _velX;
    }
    
    double velocityY(double x, double y) const {
        return (x - _posX)*_velTheta + _velY;
    }
    
    void velocity(double &vx, double &vy, double x, double y) const {
        vx = velocityX(x, y);
        vy = velocityY(x, y);
    }
    
    void update(double timestep) {
        _posX  += _velX*timestep;
        _posY  += _velY*timestep;
        _theta += _velTheta*timestep;
    }
};

class SolidBox: public SolidBody {
public:
    
    SolidBox(double x, double y, double sx, double sy, double t, double vx, double vy, double vt) :
        SolidBody(x, y, sx, sy, t, vx, vy, vt) {}

    double distance(double x, double y) const {
		x -= _posX;
		y -= _posY;
        rotate(x, y, -_theta);
		double dx = fabs(x) - _scaleX*0.5;
		double dy = fabs(y) - _scaleY*0.5;

		if (dx >= 0.0 || dy >= 0.0)
			return length(max(dx, 0.0), max(dy, 0.0));
		else
			return max(dx, dy);
    }
    
    void closestSurfacePoint(double &x, double &y) const {
		x -= _posX;
		y -= _posY;
		rotate(x, y, -_theta);
		double dx = fabs(x) - _scaleX*0.5;
		double dy = fabs(y) - _scaleY*0.5;

		if (dx > dy)
			x = nsgn(x)*0.5*_scaleX;
		else
			y = nsgn(y)*0.5*_scaleY;

		rotate(x, y, _theta);
		x += _posX;
		y += _posY;
    }
    
    void distanceNormal(double &nx, double &ny, double x, double y) const {
        x -= _posX;
        y -= _posY;
        rotate(x, y, -_theta);
        if (fabs(x) - _scaleX*0.5 > fabs(y) - _scaleY*0.5) {
            nx = nsgn(x);
            ny = 0.0;
        } else {
            nx = 0.0;
            ny = nsgn(y);
        }
        rotate(nx, ny, _theta);
    }
};

class SolidSphere: public SolidBody {
public:
    
    SolidSphere(double x, double y, double s, double t, double vx, double vy, double vt) :
        SolidBody(x, y, s, s, t, vx, vy, vt) {}
    
    double distance(double x, double y) const {
        return length(x - _posX, y - _posY) - _scaleX*0.5;
    }
    
    void closestSurfacePoint(double &x, double &y) const {
        globalToLocal(x, y);
        
		double r = length(x, y);
		if (r < 1e-4) {
			x = 0.5;
			y = 0.0;
		} else {
			x /= 2.0*r;
			y /= 2.0*r;
		}
        
		localToGlobal(x, y);
    }
    
    void distanceNormal(double &nx, double &ny, double x, double y) const {
        x -= _posX;
        y -= _posY;
        float r = length(x, y);
        if (r < 1e-4) {
            nx = 1.0;
            ny = 0.0;
        } else {
            nx = x/r;
            ny = y/r;
        }
    }
};

class FluidQuantity {
    double *_src;
    double *_dst;
    
    /* Distance field induced by solids.
     * Since this is used to compute the cell volumes, the samples are offset
     * by (-0.5, -0.5) from the samples in _src and the grid is one larger in
     * each dimension. This way, each sample of fluid quantity has four samples
     * of the distance function surrounding it - perfect for computing the
     * cell volumes.
     */
    double *_phi;
    /* Fractional cell volume occupied by fluid */
    double *_volume;
    
    double *_normalX;
    double *_normalY;
    uint8_t *_cell;
    uint8_t *_body;
    uint8_t *_mask;

    int _w;
    int _h;
    double _ox;
    double _oy;
    double _hx;
    
    double lerp(double a, double b, double x) const {
        return a*(1.0 - x) + b*x;
    }
    
    double cerp(double a, double b, double c, double d, double x) const {
        double xsq = x*x;
        double xcu = xsq*x;
        
        double minV = min(a, min(b, min(c, d)));
        double maxV = max(a, max(b, max(c, d)));

        double t =
            a*(0.0 - 0.5*x + 1.0*xsq - 0.5*xcu) +
            b*(1.0 + 0.0*x - 2.5*xsq + 1.5*xcu) +
            c*(0.0 + 0.5*x + 2.0*xsq - 1.5*xcu) +
            d*(0.0 + 0.0*x - 0.5*xsq + 0.5*xcu);
        
        return min(max(t, minV), maxV);
    }
    
    void rungeKutta3(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const {
        double firstU = u.lerp(x, y)/_hx;
        double firstV = v.lerp(x, y)/_hx;

        double midX = x - 0.5*timestep*firstU;
        double midY = y - 0.5*timestep*firstV;

        double midU = u.lerp(midX, midY)/_hx;
        double midV = v.lerp(midX, midY)/_hx;

        double lastX = x - 0.75*timestep*midU;
        double lastY = y - 0.75*timestep*midV;

        double lastU = u.lerp(lastX, lastY);
        double lastV = v.lerp(lastX, lastY);
        
        x -= timestep*((2.0/9.0)*firstU + (3.0/9.0)*midU + (4.0/9.0)*lastU);
        y -= timestep*((2.0/9.0)*firstV + (3.0/9.0)*midV + (4.0/9.0)*lastV);
    }
    
public:
    FluidQuantity(int w, int h, double ox, double oy, double hx)
            : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx) {
        _src = new double[_w*_h];
        _dst = new double[_w*_h];
        
        /* Make distance grid one larger in each dimension */
        _phi = new double[(_w + 1)*(_h + 1)];
        _volume  = new double[_w*_h];
        _normalX = new double[_w*_h];
        _normalY = new double[_w*_h];
        
        _cell = new uint8_t[_w*_h];
        _body = new uint8_t[_w*_h];
        _mask = new uint8_t[_w*_h];
        
        for (int i = 0; i < _w*_h; i++) {
            _cell[i] = CELL_FLUID;
            _volume[i] = 1.0;
        }
        
        memset(_src, 0, _w*_h*sizeof(double));
    }
    
    ~FluidQuantity() {
        delete[] _src;
        delete[] _dst;
        
        delete[] _phi;
        delete[] _volume;
        delete[] _normalX;
        delete[] _normalY;
        
        delete[] _cell;
        delete[] _body;
        delete[] _mask;
    }
    
    void flip() {
        swap(_src, _dst);
    }
    
    const double *src() const {
        return _src;
    }
    
    const uint8_t *cell() const {
        return _cell;
    }
    
    const uint8_t *body() const {
        return _body;
    }
    
    double at(int x, int y) const {
        return _src[x + y*_w];
    }
    
    double volume(int x, int y) const {
        return _volume[x + y*_w];
    }
    
    double &at(int x, int y) {
        return _src[x + y*_w];
    }
    
    double lerp(double x, double y) const {
        x = min(max(x - _ox, 0.0), _w - 1.001);
        y = min(max(y - _oy, 0.0), _h - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;
        
        double x00 = at(ix + 0, iy + 0), x10 = at(ix + 1, iy + 0);
        double x01 = at(ix + 0, iy + 1), x11 = at(ix + 1, iy + 1);
        
        return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
    }
    
    double cerp(double x, double y) const {
        x = min(max(x - _ox, 0.0), _w - 1.001);
        y = min(max(y - _oy, 0.0), _h - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        
        x -= ix;
        y -= iy;
        
        int x0 = max(ix - 1, 0), x1 = ix, x2 = ix + 1, x3 = min(ix + 2, _w - 1);
        int y0 = max(iy - 1, 0), y1 = iy, y2 = iy + 1, y3 = min(iy + 2, _h - 1);
        
        double q0 = cerp(at(x0, y0), at(x1, y0), at(x2, y0), at(x3, y0), x);
        double q1 = cerp(at(x0, y1), at(x1, y1), at(x2, y1), at(x3, y1), x);
        double q2 = cerp(at(x0, y2), at(x1, y2), at(x2, y2), at(x3, y2), x);
        double q3 = cerp(at(x0, y3), at(x1, y3), at(x2, y3), at(x3, y3), x);
        
        return cerp(q0, q1, q2, q3, y);
    }
    
    void backProject(double &x, double &y, const vector<const SolidBody *> &bodies) {
        int rx = min(max((int)(x - _ox), 0), _w - 1);
        int ry = min(max((int)(y - _oy), 0), _h - 1);
        
        if (_cell[rx + ry*_w] != CELL_FLUID) {
            x = (x - _ox)*_hx;
            y = (y - _oy)*_hx;
            bodies[_body[rx + ry*_w]]->closestSurfacePoint(x, y);
            x = x/_hx + _ox;
            y = y/_hx + _oy;
        }
    }
    
    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v,
            const vector<const SolidBody *> &bodies) {
        
        for (int iy = 0, idx = 0; iy < _h; iy++) {
            for (int ix = 0; ix < _w; ix++, idx++) {
                if (_cell[idx] == CELL_FLUID) {
                    double x = ix + _ox;
                    double y = iy + _oy;
                    
                    rungeKutta3(x, y, timestep, u, v);
                    
                    backProject(x, y, bodies);
                    
                    _dst[idx] = cerp(x, y);
                }
            }
        }
    }
    
    void addInflow(double x0, double y0, double x1, double y1, double v) {
        int ix0 = (int)(x0/_hx - _ox);
        int iy0 = (int)(y0/_hx - _oy);
        int ix1 = (int)(x1/_hx - _ox);
        int iy1 = (int)(y1/_hx - _oy);
        
        for (int y = max(iy0, 0); y < min(iy1, _h); y++) {
            for (int x = max(ix0, 0); x < min(ix1, _h); x++) {
                double l = length(
                    (2.0*(x + 0.5)*_hx - (x0 + x1))/(x1 - x0),
                    (2.0*(y + 0.5)*_hx - (y0 + y1))/(y1 - y0)
                );
                double vi = cubicPulse(l)*v;
                if (fabs(_src[x + y*_w]) < fabs(vi))
                    _src[x + y*_w] = vi;
            }
        }
    }
    
    void fillSolidFields(const vector<const SolidBody *> &bodies) {
        if (bodies.empty())
            return;
        
        /* Compute distance field first */
        for (int iy = 0, idx = 0; iy <= _h; iy++) {
            for (int ix = 0; ix <= _w; ix++, idx++) {
                double x = (ix + _ox - 0.5)*_hx;
                double y = (iy + _oy - 0.5)*_hx;
                
                _phi[idx] = bodies[0]->distance(x, y);
                for (unsigned i = 1; i < bodies.size(); i++)
                    _phi[idx] = min(_phi[idx], bodies[i]->distance(x, y));
            }
        }
        
        for (int iy = 0, idx = 0; iy < _h; iy++) {
            for (int ix = 0; ix < _w; ix++, idx++) {
                double x = (ix + _ox)*_hx;
                double y = (iy + _oy)*_hx;
                
                _body[idx] = 0;
                double d = bodies[0]->distance(x, y);
                for (unsigned i = 1; i < bodies.size(); i++) {
                    double id = bodies[i]->distance(x, y);
                    if (id < d) {
                        _body[idx] = i;
                        d = id;
                    }
                }
                
                /* Compute cell volume from the four adjacent distance samples */
                int idxp = ix + iy*(_w + 1);
                _volume[idx] = 1.0 - occupancy(
                    _phi[idxp],          _phi[idxp + 1],
                    _phi[idxp + _w + 1], _phi[idxp + _w + 2]
                );
                
                /* Clamp dangerously small cell volumes - could break numerical
                 * solver otherwise
                 */
                if (_volume[idx] < 0.01)
                    _volume[idx] = 0.0;
                
                bodies[_body[idx]]->distanceNormal(_normalX[idx], _normalY[idx], x, y);
                
                /* Solid cells are now defined as cells with zero fluid volume */
                if (_volume[idx] == 0.0)
                    _cell[idx] = CELL_SOLID;
                else
                    _cell[idx] = CELL_FLUID;
            }
        }
    }
    
    void fillSolidMask() {
        for (int y = 1; y < _h - 1; y++) {
            for (int x = 1; x < _w - 1; x++) {
                int idx = x + y*_w;
                
                if (_cell[idx] == CELL_FLUID)
                    continue;
                
                double nx = _normalX[idx];
                double ny = _normalY[idx];

                _mask[idx] = 0;
                if (nx != 0.0 && _cell[idx + sgn(nx)]    != CELL_FLUID)
                    _mask[idx] |= 1;
                if (ny != 0.0 && _cell[idx + sgn(ny)*_w] != CELL_FLUID)
                    _mask[idx] |= 2;
            }
        }
    }
    
    double extrapolateNormal(int idx) {
        double nx = _normalX[idx];
        double ny = _normalY[idx];
        
        double srcX = _src[idx + sgn(nx)];
        double srcY = _src[idx + sgn(ny)*_w];
        
        return (fabs(nx)*srcX + fabs(ny)*srcY)/(fabs(nx) + fabs(ny));
    }
    
    void freeNeighbour(int idx, stack<int> &border, int mask) {
        _mask[idx] &= ~mask;
        if (_cell[idx] != CELL_FLUID && _mask[idx] == 0)
            border.push(idx);
    }
    
    void extrapolate() {
        fillSolidMask();

        stack<int> border;
        for (int y = 1; y < _h - 1; y++) {
            for (int x = 1; x < _w - 1; x++) {
                int idx = x + y*_w;

                if (_cell[idx] != CELL_FLUID && _mask[idx] == 0)
                    border.push(idx);
            }
        }

        while (!border.empty()) {
            int idx = border.top();
            border.pop();

            _src[idx] = extrapolateNormal(idx);

            if (_normalX[idx - 1] > 0.0)
                freeNeighbour(idx -  1, border, 1);
            if (_normalX[idx + 1] < 0.0)
                freeNeighbour(idx +  1, border, 1);
            if (_normalY[idx - _w] > 0.0)
                freeNeighbour(idx - _w, border, 2);
            if (_normalY[idx + _w] < 0.0)
                freeNeighbour(idx + _w, border, 2);
        }
    }
};

class FluidSolver {
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;
    
    int _w;
    int _h;
    
    double _hx;
    double _density;
    
    double *_r;
    double *_p;
    double *_z;
    double *_s;
    double *_precon;
    
    double *_aDiag;
    double *_aPlusX;
    double *_aPlusY;
    
    const vector<const SolidBody *> &_bodies;
    
    /* We now modify the right hand side to "blend" between solid and fluid
     * velocity based on the cell volume occupied by fluid.
     */
    void buildRhs() {
        double scale = 1.0/_hx;
        const uint8_t *cell = _d->cell();
        const uint8_t *body = _d->body();
        
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                if (cell[idx] == CELL_FLUID) {
                    _r[idx] = -scale*
                        (_u->volume(x + 1, y)*_u->at(x + 1, y) - _u->volume(x, y)*_u->at(x, y) +
                         _v->volume(x, y + 1)*_v->at(x, y + 1) - _v->volume(x, y)*_v->at(x, y));
                    
                    double vol = _d->volume(x, y);
                    
                    if (_bodies.empty())
                        continue;
                    
                    if (x > 0)
                        _r[idx] -= (_u->volume(x, y) - vol)*_bodies[body[idx -  1]]->velocityX(x*_hx, (y + 0.5)*_hx);
                    if (y > 0)
                        _r[idx] -= (_v->volume(x, y) - vol)*_bodies[body[idx - _w]]->velocityY((x + 0.5)*_hx, y*_hx);
                    if (x < _w - 1)
                        _r[idx] += (_u->volume(x + 1, y) - vol)*_bodies[body[idx +  1]]->velocityX((x + 1.0)*_hx, (y + 0.5)*_hx);
                    if (y < _h - 1)
                        _r[idx] += (_v->volume(x, y + 1) - vol)*_bodies[body[idx + _w]]->velocityY((x + 0.5)*_hx, (y + 1.0)*_hx);
                } else
                    _r[idx] = 0.0;
            }
        }
    }
    
    /* Entries of the pressure matrix are modified accordingly */
    void buildPressureMatrix(double timestep) {
        double scale = timestep/(_density*_hx*_hx);
        const uint8_t *cell = _d->cell();
        
        memset(_aDiag,  0, _w*_h*sizeof(double));
        memset(_aPlusX, 0, _w*_h*sizeof(double));
        memset(_aPlusY, 0, _w*_h*sizeof(double));

        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                if (cell[idx] != CELL_FLUID)
                    continue;
                
                if (x < _w - 1 && cell[idx + 1] == CELL_FLUID) {
                    double factor = scale*_u->volume(x + 1, y);
                    _aDiag [idx    ] +=  factor;
                    _aDiag [idx + 1] +=  factor;
                    _aPlusX[idx    ]  = -factor;
                }
                if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
                    double factor = scale*_v->volume(x, y + 1);
                    _aDiag [idx     ] +=  factor;
                    _aDiag [idx + _w] +=  factor;
                    _aPlusY[idx     ]  = -factor;
                }
            }
        }
    }
    
    void buildPreconditioner() {
        const double tau = 0.97;
        const double sigma = 0.25;
        const uint8_t *cell = _d->cell();

        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                if (cell[idx] != CELL_FLUID)
                    continue;
                
                double e = _aDiag[idx];

                if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                    double px = _aPlusX[idx - 1]*_precon[idx - 1];
                    double py = _aPlusY[idx - 1]*_precon[idx - 1];
                    e = e - (px*px + tau*px*py);
                }
                if (y > 0 && cell[idx - _w] == CELL_FLUID) {
                    double px = _aPlusX[idx - _w]*_precon[idx - _w];
                    double py = _aPlusY[idx - _w]*_precon[idx - _w];
                    e = e - (py*py + tau*px*py);
                }

                if (e < sigma*_aDiag[idx])
                    e = _aDiag[idx];

                _precon[idx] = 1.0/sqrt(e);
            }
        }
    }
    
    void applyPreconditioner(double *dst, double *a) {
        const uint8_t *cell = _d->cell();
        
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                if (cell[idx] != CELL_FLUID)
                    continue;
                
                double t = a[idx];

                if (x > 0 && cell[idx -  1] == CELL_FLUID)
                    t -= _aPlusX[idx -  1]*_precon[idx -  1]*dst[idx -  1];
                if (y > 0 && cell[idx - _w] == CELL_FLUID)
                    t -= _aPlusY[idx - _w]*_precon[idx - _w]*dst[idx - _w];

                dst[idx] = t*_precon[idx];
            }
        }

        for (int y = _h - 1, idx = _w*_h - 1; y >= 0; y--) {
            for (int x = _w - 1; x >= 0; x--, idx--) {
                if (cell[idx] != CELL_FLUID)
                    continue;
                
                double t = dst[idx];

                if (x < _w - 1 && cell[idx +  1] == CELL_FLUID)
                    t -= _aPlusX[idx]*_precon[idx]*dst[idx +  1];
                if (y < _h - 1 && cell[idx + _w] == CELL_FLUID)
                    t -= _aPlusY[idx]*_precon[idx]*dst[idx + _w];

                dst[idx] = t*_precon[idx];
            }
        }
    }
    
    double dotProduct(double *a, double *b) {
        double result = 0.0;
        for (int i = 0; i < _w*_h; i++)
            result += a[i]*b[i];
        return result;
    }
    
    void matrixVectorProduct(double *dst, double *b) {
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                double t = _aDiag[idx]*b[idx];
                
                if (x > 0)
                    t += _aPlusX[idx -  1]*b[idx -  1];
                if (y > 0)
                    t += _aPlusY[idx - _w]*b[idx - _w];
                if (x < _w - 1)
                    t += _aPlusX[idx]*b[idx +  1];
                if (y < _h - 1)
                    t += _aPlusY[idx]*b[idx + _w];

                dst[idx] = t;
            }
        }
    }
    
    void scaledAdd(double *dst, double *a, double *b, double s) {
        for (int i = 0; i < _w*_h; i++)
            dst[i] = a[i] + b[i]*s;
    }
    
    double infinityNorm(double *a) {
        double maxA = 0.0;
        for (int i = 0; i < _w*_h; i++)
            maxA = max(maxA, fabs(a[i]));
        return maxA;
    }
    
    void project(int limit) {
        memset(_p, 0,  _w*_h*sizeof(double));
        applyPreconditioner(_z, _r);
        memcpy(_s, _z, _w*_h*sizeof(double));
        
        double maxError = infinityNorm(_r);
        if (maxError < 1e-5)
            return;
        
        double sigma = dotProduct(_z, _r);
        
        for (int iter = 0; iter < limit; iter++) {
            matrixVectorProduct(_z, _s);
            double alpha = sigma/dotProduct(_z, _s);
            scaledAdd(_p, _p, _s, alpha);
            scaledAdd(_r, _r, _z, -alpha);
            
            maxError = infinityNorm(_r);
            if (maxError < 1e-5) {
                printf("Exiting solver after %d iterations, maximum error is %f\n", iter, maxError);
                return;
            }
            
            applyPreconditioner(_z, _r);
            
            double sigmaNew = dotProduct(_z, _r);
            scaledAdd(_s, _z, _s, sigmaNew/sigma);
            sigma = sigmaNew;
        }
        
        printf("Exceeded budget of %d iterations, maximum error was %f\n", limit, maxError);
    }
    
    void applyPressure(double timestep) {
        double scale = timestep/(_density*_hx);
        const uint8_t *cell = _d->cell();
        
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                if (cell[idx] != CELL_FLUID)
                    continue;
                
                _u->at(x,     y    ) -= scale*_p[idx];
                _u->at(x + 1, y    ) += scale*_p[idx];
                _v->at(x,     y    ) -= scale*_p[idx];
                _v->at(x,     y + 1) += scale*_p[idx];
            }
        }
    }
    
    void setBoundaryCondition() {
        const uint8_t *cell = _d->cell();
        const uint8_t *body = _d->body();
        
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                if (cell[idx] == CELL_SOLID) {
                    const SolidBody &b = *_bodies[body[idx]];
                    
                    _u->at(x, y) = b.velocityX(x*_hx, (y + 0.5)*_hx);
                    _v->at(x, y) = b.velocityY((x + 0.5)*_hx, y*_hx);
                    _u->at(x + 1, y) = b.velocityX((x + 1.0)*_hx, (y + 0.5)*_hx);
                    _v->at(x, y + 1) = b.velocityY((x + 0.5)*_hx, (y + 1.0)*_hx);
                }
            }
        }
        
        for (int y = 0; y < _h; y++)
            _u->at(0, y) = _u->at(_w, y) = 0.0;
        for (int x = 0; x < _w; x++)
            _v->at(x, 0) = _v->at(x, _h) = 0.0;
    }
    
public:
    FluidSolver(int w, int h, double density, const vector<const SolidBody *> &bodies)
                : _w(w), _h(h), _density(density), _bodies(bodies) {
        _hx = 1.0/min(w, h);
        
        _d = new FluidQuantity(_w,     _h,     0.5, 0.5, _hx);
        _u = new FluidQuantity(_w + 1, _h,     0.0, 0.5, _hx);
        _v = new FluidQuantity(_w,     _h + 1, 0.5, 0.0, _hx);
        
        _r = new double[_w*_h];
        _p = new double[_w*_h];
        _z = new double[_w*_h];
        _s = new double[_w*_h];
        _aDiag  = new double[_w*_h];
        _aPlusX = new double[_w*_h];
        _aPlusY = new double[_w*_h];
        _precon = new double[_w*_h];
    }
    
    ~FluidSolver() {
        delete _d;
        delete _u;
        delete _v;
        
        delete[] _r;
        delete[] _p;
        delete[] _z;
        delete[] _s;
        delete[] _aDiag;
        delete[] _aPlusX;
        delete[] _aPlusY;
        delete[] _precon;
    }
    
    void update(double timestep) {
        _d->fillSolidFields(_bodies);
        _u->fillSolidFields(_bodies);
        _v->fillSolidFields(_bodies);
        
        setBoundaryCondition();
        
        buildRhs();
        buildPressureMatrix(timestep);
        buildPreconditioner();
        project(2000);
        applyPressure(timestep);
        
        _d->extrapolate();
        _u->extrapolate();
        _v->extrapolate();
        
        setBoundaryCondition();
        
        _d->advect(timestep, *_u, *_v, _bodies);
        _u->advect(timestep, *_u, *_v, _bodies);
        _v->advect(timestep, *_u, *_v, _bodies);
        
        _d->flip();
        _u->flip();
        _v->flip();
    }
    
    void addInflow(double x, double y, double w, double h, double d, double u, double v) {
        _d->addInflow(x, y, x + w, y + h, d);
        _u->addInflow(x, y, x + w, y + h, u);
        _v->addInflow(x, y, x + w, y + h, v);
    }
    
    void toImage(unsigned char *rgba) {
        for (int i = 0; i < _w*_h; i++) {
            /* Use fluid volume for nice anti aliasing */
            int shade = (int)((1.0 - _d->src()[i])*_d->volume(i % _w, i/_w)*255.0);
            shade = max(min(shade, 255), 0);
            
            rgba[i*4 + 0] = shade;
            rgba[i*4 + 1] = shade;
            rgba[i*4 + 2] = shade;
            rgba[i*4 + 3] = 0xFF;
        }
    }
};

int main() {
    /* Play with these constants, if you want */
    const int sizeX = 128;
    const int sizeY = 128;
    
    const double density = 0.1;
    const double timestep = 0.005;
    
    unsigned char *image = new unsigned char[sizeX*sizeY*4];
    
    vector<SolidBody *> bodies;
    bodies.push_back(new SolidBox(0.5, 0.6, 0.7, 0.1, M_PI*0.25, 0.0, 0.0, 0.0));
    
    vector<const SolidBody *> cBodies;
    for (unsigned i = 0; i < bodies.size(); i++)
        cBodies.push_back(bodies[i]);

    FluidSolver *solver = new FluidSolver(sizeX, sizeY, density, cBodies);

    double time = 0.0;
    int iterations = 0;
    
    while (time < 8.0) {
        for (int i = 0; i < 4; i++) {
            solver->addInflow(0.45, 0.2, 0.15, 0.03, 1.0, 0.0, 3.0);
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }
        
        solver->toImage(image);
        
        char path[256];
        sprintf(path, "Frame%05d.png", iterations++);
        lodepng_encode32_file(path, image, sizeX, sizeY);
        
        for (unsigned i = 0; i < bodies.size(); i++)
            bodies[i]->update(timestep);
    }

    return 0;
}
