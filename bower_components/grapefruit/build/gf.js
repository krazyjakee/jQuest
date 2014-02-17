/**
 * @license
 * GrapeFruit Game Engine - v0.2.0
 * Copyright © 2012-2014, Chad Engler
 * https://github.com/grapefruitjs/grapefruit
 *
 * Compiled: 2014-02-17
 *
 * GrapeFruit Game Engine is licensed under the MIT License.
 * http://www.opensource.org/licenses/mit-license.php
 */
!function(e){if("object"==typeof exports)module.exports=e();else if("function"==typeof define&&define.amd)define(e);else{var f;"undefined"!=typeof window?f=window:"undefined"!=typeof global?f=global:"undefined"!=typeof self&&(f=self),f.gf=e()}}(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);throw new Error("Cannot find module '"+o+"'")}var f=n[o]={exports:{}};t[o][0].call(f.exports,function(e){var n=t[o][1][e];return s(n?n:e)},f,f.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(_dereq_,module,exports){
(function(){
/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

Object.create = Object.create || function(o) {
	function F() {}
	F.prototype = o;
	return new F();
};
 
//var VERSION = CP_VERSION_MAJOR + "." + CP_VERSION_MINOR + "." + CP_VERSION_RELEASE;

var cp;
if(typeof exports === 'undefined'){
	cp = {};

	if(typeof window === 'object'){
		window.cp = cp;
	}
} else {
	cp = exports;
}

var assert = function(value, message)
{
	if (!value) {
		throw new Error('Assertion failed: ' + message);
	}
};

var assertSoft = function(value, message)
{
	if(!value && console && console.warn) {
		console.warn("ASSERTION FAILED: " + message);
		if(console.trace) {
			console.trace();
		}
	}
};

var mymin = function(a, b)
{
	return a < b ? a : b;
};
var mymax = function(a, b)
{
	return a > b ? a : b;
};

var min, max;
if (typeof window === 'object' && window.navigator.userAgent.indexOf('Firefox') > -1){
	// On firefox, Math.min and Math.max are really fast:
	// http://jsperf.com/math-vs-greater-than/8
	min = Math.min;
	max = Math.max;
} else {
	// On chrome and safari, Math.min / max are slooow. The ternery operator above is faster
	// than the builtins because we only have to deal with 2 arguments that are always numbers.
	min = mymin;
	max = mymax;
}

/* The hashpair function takes two numbers and returns a hash code for them.
 * Required that hashPair(a, b) === hashPair(b, a).
 * Chipmunk's hashPair function is defined as:
 *   #define CP_HASH_COEF (3344921057ul)
 *   #define CP_HASH_PAIR(A, B) ((cpHashValue)(A)*CP_HASH_COEF ^ (cpHashValue)(B)*CP_HASH_COEF)
 * But thats not suitable in javascript because multiplying by a large number will make the number
 * a large float.
 *
 * The result of hashPair is used as the key in objects, so it returns a string.
 */
var hashPair = function(a, b)
{
	//assert(typeof(a) === 'number', "HashPair used on something not a number");
	return a < b ? a + ' ' + b : b + ' ' + a;
};

var deleteObjFromList = function(arr, obj)
{
	for(var i=0; i<arr.length; i++){
		if(arr[i] === obj){
			arr[i] = arr[arr.length - 1];
			arr.length--;
			
			return;
		}
	}
};

var closestPointOnSegment = function(p, a, b)
{
	var delta = vsub(a, b);
	var t = clamp01(vdot(delta, vsub(p, b))/vlengthsq(delta));
	return vadd(b, vmult(delta, t));
};

var closestPointOnSegment2 = function(px, py, ax, ay, bx, by)
{
	var deltax = ax - bx;
	var deltay = ay - by;
	var t = clamp01(vdot2(deltax, deltay, px - bx, py - by)/vlengthsq2(deltax, deltay));
	return new Vect(bx + deltax * t, by + deltay * t);
};

var momentForCircle = cp.momentForCircle = function(m, r1, r2, offset)
{
	return m*(0.5*(r1*r1 + r2*r2) + vlengthsq(offset));
};

var areaForCircle = cp.areaForCircle = function(r1, r2)
{
	return Math.PI*Math.abs(r1*r1 - r2*r2);
};

var momentForSegment = cp.momentForSegment = function(m, a, b)
{
	var offset = vmult(vadd(a, b), 0.5);
	return m*(vdistsq(b, a)/12 + vlengthsq(offset));
};

var areaForSegment = cp.areaForSegment = function(a, b, r)
{
	return r*(Math.PI*r + 2*vdist(a, b));
};

var momentForPoly = cp.momentForPoly = function(m, verts, offset)
{
	var sum1 = 0;
	var sum2 = 0;
	var len = verts.length;
	for(var i=0; i<len; i+=2){
		var v1x = verts[i] + offset.x;
	 	var v1y = verts[i+1] + offset.y;
		var v2x = verts[(i+2)%len] + offset.x;
		var v2y = verts[(i+3)%len] + offset.y;

		var a = vcross2(v2x, v2y, v1x, v1y);
		var b = vdot2(v1x, v1y, v1x, v1y) + vdot2(v1x, v1y, v2x, v2y) + vdot2(v2x, v2y, v2x, v2y);
		
		sum1 += a*b;
		sum2 += a;
	}
	
	return (m*sum1)/(6*sum2);
};

var areaForPoly = cp.areaForPoly = function(verts)
{
	var area = 0;
	for(var i=0, len=verts.length; i<len; i+=2){
		area += vcross(new Vect(verts[i], verts[i+1]), new Vect(verts[(i+2)%len], verts[(i+3)%len]));
	}
	
	return -area/2;
};

var centroidForPoly = cp.centroidForPoly = function(verts)
{
	var sum = 0;
	var vsum = new Vect(0,0);
	
	for(var i=0, len=verts.length; i<len; i+=2){
		var v1 = new Vect(verts[i], verts[i+1]);
		var v2 = new Vect(verts[(i+2)%len], verts[(i+3)%len]);
		var cross = vcross(v1, v2);
		
		sum += cross;
		vsum = vadd(vsum, vmult(vadd(v1, v2), cross));
	}
	
	return vmult(vsum, 1/(3*sum));
};

var recenterPoly = cp.recenterPoly = function(verts)
{
	var centroid = centroidForPoly(verts);
	
	for(var i=0; i<verts.length; i+=2){
		verts[i] -= centroid.x;
		verts[i+1] -= centroid.y;
	}
};

var momentForBox = cp.momentForBox = function(m, width, height)
{
	return m*(width*width + height*height)/12;
};

var momentForBox2 = cp.momentForBox2 = function(m, box)
{
	width = box.r - box.l;
	height = box.t - box.b;
	offset = vmult([box.l + box.r, box.b + box.t], 0.5);
	
	// TODO NaN when offset is 0 and m is INFINITY	
	return momentForBox(m, width, height) + m*vlengthsq(offset);
};

// Quick hull

var loopIndexes = cp.loopIndexes = function(verts)
{
	var start = 0, end = 0;
	var minx, miny, maxx, maxy;
	minx = maxx = verts[0];
	miny = maxy = verts[1];
	
	var count = verts.length >> 1;
  for(var i=1; i<count; i++){
		var x = verts[i*2];
		var y = verts[i*2 + 1];
		
    if(x < minx || (x == minx && y < miny)){
			minx = x;
			miny = y;
      start = i;
    } else if(x > maxx || (x == maxx && y > maxy)){
			maxx = x;
			maxy = y;
			end = i;
		}
	}
	return [start, end];
};

var SWAP = function(arr, idx1, idx2)
{
	var tmp = arr[idx1*2];
	arr[idx1*2] = arr[idx2*2];
	arr[idx2*2] = tmp;

	tmp = arr[idx1*2+1];
	arr[idx1*2+1] = arr[idx2*2+1];
	arr[idx2*2+1] = tmp;
};

var QHullPartition = function(verts, offs, count, a, b, tol)
{
	if(count === 0) return 0;
	
	var max = 0;
	var pivot = offs;
	
	var delta = vsub(b, a);
	var valueTol = tol * vlength(delta);
	
	var head = offs;
	for(var tail = offs+count-1; head <= tail;){
		var v = new Vect(verts[head * 2], verts[head * 2 + 1]);
		var value = vcross(delta, vsub(v, a));
		if(value > valueTol){
			if(value > max){
				max = value;
				pivot = head;
			}
			
			head++;
		} else {
			SWAP(verts, head, tail);
			tail--;
		}
	}
	
	// move the new pivot to the front if it's not already there.
	if(pivot != offs) SWAP(verts, offs, pivot);
	return head - offs;
};

var QHullReduce = function(tol, verts, offs, count, a, pivot, b, resultPos)
{
	if(count < 0){
		return 0;
	} else if(count == 0) {
		verts[resultPos*2] = pivot.x;
		verts[resultPos*2+1] = pivot.y;
		return 1;
	} else {
		var left_count = QHullPartition(verts, offs, count, a, pivot, tol);
		var left = new Vect(verts[offs*2], verts[offs*2+1]);
		var index = QHullReduce(tol, verts, offs + 1, left_count - 1, a, left, pivot, resultPos);
		
		var pivotPos = resultPos + index++;
		verts[pivotPos*2] = pivot.x;
		verts[pivotPos*2+1] = pivot.y;
		
		var right_count = QHullPartition(verts, offs + left_count, count - left_count, pivot, b, tol);
		var right = new Vect(verts[(offs+left_count)*2], verts[(offs+left_count)*2+1]);
		return index + QHullReduce(tol, verts, offs + left_count + 1, right_count - 1, pivot, right, b, resultPos + index);
	}
};

// QuickHull seemed like a neat algorithm, and efficient-ish for large input sets.
// My implementation performs an in place reduction using the result array as scratch space.
//
// Pass an Array into result to put the result of the calculation there. Otherwise, pass null
// and the verts list will be edited in-place.
//
// Expects the verts to be described in the same way as cpPolyShape - which is to say, it should
// be a list of [x1,y1,x2,y2,x3,y3,...].
var convexHull = cp.convexHull = function(verts, result, tol)
{
	if(result){
		// Copy the line vertexes into the empty part of the result polyline to use as a scratch buffer.
		for (var i = 0; i < verts.length; i++){
			result[i] = verts[i];
		}
	} else {
		// If a result array was not specified, reduce the input instead.
		result = verts;
	}
	
	// Degenerate case, all points are the same.
	var indexes = loopIndexes(verts);
	var start = indexes[0], end = indexes[1];
	if(start == end){
		//if(first) (*first) = 0;
		result.length = 2;
		return result;
	}
	
	SWAP(result, 0, start);
	SWAP(result, 1, end == 0 ? start : end);
	
	var a = new Vect(result[0], result[1]);
	var b = new Vect(result[2], result[3]);
	
	var count = verts.length >> 1;
	//if(first) (*first) = start;
	var resultCount = QHullReduce(tol, result, 2, count - 2, a, b, a, 1) + 1;
	result.length = resultCount*2;

	assertSoft(polyValidate(result),
		"Internal error: cpConvexHull() and cpPolyValidate() did not agree." +
		"Please report this error with as much info as you can.");
	return result;
};

/// Clamp @c f to be between @c min and @c max.
var clamp = function(f, minv, maxv)
{
	return min(max(f, minv), maxv);
};

/// Clamp @c f to be between 0 and 1.
var clamp01 = function(f)
{
	return max(0, min(f, 1));
};

/// Linearly interpolate (or extrapolate) between @c f1 and @c f2 by @c t percent.
var lerp = function(f1, f2, t)
{
	return f1*(1 - t) + f2*t;
};

/// Linearly interpolate from @c f1 to @c f2 by no more than @c d.
var lerpconst = function(f1, f2, d)
{
	return f1 + clamp(f2 - f1, -d, d);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// I'm using an array tuple here because (at time of writing) its about 3x faster
// than an object on firefox, and the same speed on chrome.

var numVects = 0;

var traces = {};

var Vect = cp.Vect = function(x, y)
{
	this.x = x;
	this.y = y;
	numVects++;

//	var s = new Error().stack;
//	traces[s] = traces[s] ? traces[s]+1 : 1;
};

cp.v = function (x,y) { return new Vect(x, y) };

var vzero = cp.vzero = new Vect(0,0);

// The functions below *could* be rewritten to be instance methods on Vect. I don't
// know how that would effect performance. For now, I'm keeping the JS similar to
// the original C code.

/// Vector dot product.
var vdot = cp.v.dot = function(v1, v2)
{
	return v1.x*v2.x + v1.y*v2.y;
};

var vdot2 = function(x1, y1, x2, y2)
{
	return x1*x2 + y1*y2;
};

/// Returns the length of v.
var vlength = cp.v.len = function(v)
{
	return Math.sqrt(vdot(v, v));
};

var vlength2 = cp.v.len2 = function(x, y)
{
	return Math.sqrt(x*x + y*y);
};

/// Check if two vectors are equal. (Be careful when comparing floating point numbers!)
var veql = cp.v.eql = function(v1, v2)
{
	return (v1.x === v2.x && v1.y === v2.y);
};

/// Add two vectors
var vadd = cp.v.add = function(v1, v2)
{
	return new Vect(v1.x + v2.x, v1.y + v2.y);
};

Vect.prototype.add = function(v2)
{
	this.x += v2.x;
	this.y += v2.y;
	return this;
};

/// Subtract two vectors.
var vsub = cp.v.sub = function(v1, v2)
{
	return new Vect(v1.x - v2.x, v1.y - v2.y);
};

Vect.prototype.sub = function(v2)
{
	this.x -= v2.x;
	this.y -= v2.y;
	return this;
};

/// Negate a vector.
var vneg = cp.v.neg = function(v)
{
	return new Vect(-v.x, -v.y);
};

Vect.prototype.neg = function()
{
	this.x = -this.x;
	this.y = -this.y;
	return this;
};

/// Scalar multiplication.
var vmult = cp.v.mult = function(v, s)
{
	return new Vect(v.x*s, v.y*s);
};

Vect.prototype.mult = function(s)
{
	this.x *= s;
	this.y *= s;
	return this;
};

/// 2D vector cross product analog.
/// The cross product of 2D vectors results in a 3D vector with only a z component.
/// This function returns the magnitude of the z value.
var vcross = cp.v.cross = function(v1, v2)
{
	return v1.x*v2.y - v1.y*v2.x;
};

var vcross2 = function(x1, y1, x2, y2)
{
	return x1*y2 - y1*x2;
};

/// Returns a perpendicular vector. (90 degree rotation)
var vperp = cp.v.perp = function(v)
{
	return new Vect(-v.y, v.x);
};

/// Returns a perpendicular vector. (-90 degree rotation)
var vpvrperp = cp.v.pvrperp = function(v)
{
	return new Vect(v.y, -v.x);
};

/// Returns the vector projection of v1 onto v2.
var vproject = cp.v.project = function(v1, v2)
{
	return vmult(v2, vdot(v1, v2)/vlengthsq(v2));
};

Vect.prototype.project = function(v2)
{
	this.mult(vdot(this, v2) / vlengthsq(v2));
	return this;
};

/// Uses complex number multiplication to rotate v1 by v2. Scaling will occur if v1 is not a unit vector.
var vrotate = cp.v.rotate = function(v1, v2)
{
	return new Vect(v1.x*v2.x - v1.y*v2.y, v1.x*v2.y + v1.y*v2.x);
};

Vect.prototype.rotate = function(v2)
{
	this.x = this.x * v2.x - this.y * v2.y;
	this.y = this.x * v2.y + this.y * v2.x;
	return this;
};

/// Inverse of vrotate().
var vunrotate = cp.v.unrotate = function(v1, v2)
{
	return new Vect(v1.x*v2.x + v1.y*v2.y, v1.y*v2.x - v1.x*v2.y);
};

/// Returns the squared length of v. Faster than vlength() when you only need to compare lengths.
var vlengthsq = cp.v.lengthsq = function(v)
{
	return vdot(v, v);
};

var vlengthsq2 = cp.v.lengthsq2 = function(x, y)
{
	return x*x + y*y;
};

/// Linearly interpolate between v1 and v2.
var vlerp = cp.v.lerp = function(v1, v2, t)
{
	return vadd(vmult(v1, 1 - t), vmult(v2, t));
};

/// Returns a normalized copy of v.
var vnormalize = cp.v.normalize = function(v)
{
	return vmult(v, 1/vlength(v));
};

/// Returns a normalized copy of v or vzero if v was already vzero. Protects against divide by zero errors.
var vnormalize_safe = cp.v.normalize_safe = function(v)
{
	return (v.x === 0 && v.y === 0 ? vzero : vnormalize(v));
};

/// Clamp v to length len.
var vclamp = cp.v.clamp = function(v, len)
{
	return (vdot(v,v) > len*len) ? vmult(vnormalize(v), len) : v;
};

/// Linearly interpolate between v1 towards v2 by distance d.
var vlerpconst = cp.v.lerpconst = function(v1, v2, d)
{
	return vadd(v1, vclamp(vsub(v2, v1), d));
};

/// Returns the distance between v1 and v2.
var vdist = cp.v.dist = function(v1, v2)
{
	return vlength(vsub(v1, v2));
};

/// Returns the squared distance between v1 and v2. Faster than vdist() when you only need to compare distances.
var vdistsq = cp.v.distsq = function(v1, v2)
{
	return vlengthsq(vsub(v1, v2));
};

/// Returns true if the distance between v1 and v2 is less than dist.
var vnear = cp.v.near = function(v1, v2, dist)
{
	return vdistsq(v1, v2) < dist*dist;
};

/// Spherical linearly interpolate between v1 and v2.
var vslerp = cp.v.slerp = function(v1, v2, t)
{
	var omega = Math.acos(vdot(v1, v2));
	
	if(omega) {
		var denom = 1/Math.sin(omega);
		return vadd(vmult(v1, Math.sin((1 - t)*omega)*denom), vmult(v2, Math.sin(t*omega)*denom));
	} else {
		return v1;
	}
};

/// Spherical linearly interpolate between v1 towards v2 by no more than angle a radians
var vslerpconst = cp.v.slerpconst = function(v1, v2, a)
{
	var angle = Math.acos(vdot(v1, v2));
	return vslerp(v1, v2, min(a, angle)/angle);
};

/// Returns the unit length vector for the given angle (in radians).
var vforangle = cp.v.forangle = function(a)
{
	return new Vect(Math.cos(a), Math.sin(a));
};

/// Returns the angular direction v is pointing in (in radians).
var vtoangle = cp.v.toangle = function(v)
{
	return Math.atan2(v.y, v.x);
};

///	Returns a string representation of v. Intended mostly for debugging purposes and not production use.
var vstr = cp.v.str = function(v)
{
	return "(" + v.x.toFixed(3) + ", " + v.y.toFixed(3) + ")";
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/// Chipmunk's axis-aligned 2D bounding box type along with a few handy routines.

var numBB = 0;

// Bounding boxes are JS objects with {l, b, r, t} = left, bottom, right, top, respectively.
var BB = cp.BB = function(l, b, r, t)
{
	this.l = l;
	this.b = b;
	this.r = r;
	this.t = t;

	numBB++;
};

cp.bb = function(l, b, r, t) { return new BB(l, b, r, t); };

var bbNewForCircle = function(p, r)
{
	return new BB(
			p.x - r,
			p.y - r,
			p.x + r,
			p.y + r
		);
};

/// Returns true if @c a and @c b intersect.
var bbIntersects = function(a, b)
{
	return (a.l <= b.r && b.l <= a.r && a.b <= b.t && b.b <= a.t);
};
var bbIntersects2 = function(bb, l, b, r, t)
{
	return (bb.l <= r && l <= bb.r && bb.b <= t && b <= bb.t);
};

/// Returns true if @c other lies completely within @c bb.
var bbContainsBB = function(bb, other)
{
	return (bb.l <= other.l && bb.r >= other.r && bb.b <= other.b && bb.t >= other.t);
};

/// Returns true if @c bb contains @c v.
var bbContainsVect = function(bb, v)
{
	return (bb.l <= v.x && bb.r >= v.x && bb.b <= v.y && bb.t >= v.y);
};
var bbContainsVect2 = function(l, b, r, t, v)
{
	return (l <= v.x && r >= v.x && b <= v.y && t >= v.y);
};

/// Returns a bounding box that holds both bounding boxes.
var bbMerge = function(a, b){
	return new BB(
			min(a.l, b.l),
			min(a.b, b.b),
			max(a.r, b.r),
			max(a.t, b.t)
		);
};

/// Returns a bounding box that holds both @c bb and @c v.
var bbExpand = function(bb, v){
	return new BB(
			min(bb.l, v.x),
			min(bb.b, v.y),
			max(bb.r, v.x),
			max(bb.t, v.y)
		);
};

/// Returns the area of the bounding box.
var bbArea = function(bb)
{
	return (bb.r - bb.l)*(bb.t - bb.b);
};

/// Merges @c a and @c b and returns the area of the merged bounding box.
var bbMergedArea = function(a, b)
{
	return (max(a.r, b.r) - min(a.l, b.l))*(max(a.t, b.t) - min(a.b, b.b));
};

var bbMergedArea2 = function(bb, l, b, r, t)
{
	return (max(bb.r, r) - min(bb.l, l))*(max(bb.t, t) - min(bb.b, b));
};

/// Return true if the bounding box intersects the line segment with ends @c a and @c b.
var bbIntersectsSegment = function(bb, a, b)
{
	return (bbSegmentQuery(bb, a, b) != Infinity);
};

/// Clamp a vector to a bounding box.
var bbClampVect = function(bb, v)
{
	var x = min(max(bb.l, v.x), bb.r);
	var y = min(max(bb.b, v.y), bb.t);
	return new Vect(x, y);
};

// TODO edge case issue
/// Wrap a vector to a bounding box.
var bbWrapVect = function(bb, v)
{
	var ix = Math.abs(bb.r - bb.l);
	var modx = (v.x - bb.l) % ix;
	var x = (modx > 0) ? modx : modx + ix;
	
	var iy = Math.abs(bb.t - bb.b);
	var mody = (v.y - bb.b) % iy;
	var y = (mody > 0) ? mody : mody + iy;
	
	return new Vect(x + bb.l, y + bb.b);
};
/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
 
/// Segment query info struct.
/* These are created using literals where needed.
typedef struct cpSegmentQueryInfo {
	/// The shape that was hit, null if no collision occured.
	cpShape *shape;
	/// The normalized distance along the query segment in the range [0, 1].
	cpFloat t;
	/// The normal of the surface hit.
	cpVect n;
} cpSegmentQueryInfo;
*/

var shapeIDCounter = 0;

var CP_NO_GROUP = cp.NO_GROUP = 0;
var CP_ALL_LAYERS = cp.ALL_LAYERS = ~0;

cp.resetShapeIdCounter = function()
{
	shapeIDCounter = 0;
};

/// The cpShape struct defines the shape of a rigid body.
//
/// Opaque collision shape struct. Do not create directly - instead use
/// PolyShape, CircleShape and SegmentShape.
var Shape = cp.Shape = function(body) {
	/// The rigid body this collision shape is attached to.
	this.body = body;

	/// The current bounding box of the shape.
	this.bb_l = this.bb_b = this.bb_r = this.bb_t = 0;

	this.hashid = shapeIDCounter++;

	/// Sensor flag.
	/// Sensor shapes call collision callbacks but don't produce collisions.
	this.sensor = false;
	
	/// Coefficient of restitution. (elasticity)
	this.e = 0;
	/// Coefficient of friction.
	this.u = 0;
	/// Surface velocity used when solving for friction.
	this.surface_v = vzero;
	
	/// Collision type of this shape used when picking collision handlers.
	this.collision_type = 0;
	/// Group of this shape. Shapes in the same group don't collide.
	this.group = 0;
	// Layer bitmask for this shape. Shapes only collide if the bitwise and of their layers is non-zero.
	this.layers = CP_ALL_LAYERS;
	
	this.space = null;

	// Copy the collision code from the prototype into the actual object. This makes collision
	// function lookups slightly faster.
	this.collisionCode = this.collisionCode;
};

Shape.prototype.setElasticity = function(e) { this.e = e; };
Shape.prototype.setFriction = function(u) { this.body.activate(); this.u = u; };
Shape.prototype.setLayers = function(layers) { this.body.activate(); this.layers = layers; };
Shape.prototype.setSensor = function(sensor) { this.body.activate(); this.sensor = sensor; };
Shape.prototype.setCollisionType = function(collision_type) { this.body.activate(); this.collision_type = collision_type; };
Shape.prototype.getBody = function() { return this.body; };

Shape.prototype.active = function()
{
// return shape->prev || (shape->body && shape->body->shapeList == shape);
	return this.body && this.body.shapeList.indexOf(this) !== -1;
};

Shape.prototype.setBody = function(body)
{
	assert(!this.active(), "You cannot change the body on an active shape. You must remove the shape from the space before changing the body.");
	this.body = body;
};

Shape.prototype.cacheBB = function()
{
	return this.update(this.body.p, this.body.rot);
};

Shape.prototype.update = function(pos, rot)
{
	assert(!isNaN(rot.x), 'Rotation is NaN');
	assert(!isNaN(pos.x), 'Position is NaN');
	this.cacheData(pos, rot);
};

Shape.prototype.pointQuery = function(p)
{
	var info = this.nearestPointQuery(p);
	if (info.d < 0) return info;
};

Shape.prototype.getBB = function()
{
	return new BB(this.bb_l, this.bb_b, this.bb_r, this.bb_t);
};

/* Not implemented - all these getters and setters. Just edit the object directly.
CP_DefineShapeStructGetter(cpBody*, body, Body);
void cpShapeSetBody(cpShape *shape, cpBody *body);

CP_DefineShapeStructGetter(cpBB, bb, BB);
CP_DefineShapeStructProperty(cpBool, sensor, Sensor, cpTrue);
CP_DefineShapeStructProperty(cpFloat, e, Elasticity, cpFalse);
CP_DefineShapeStructProperty(cpFloat, u, Friction, cpTrue);
CP_DefineShapeStructProperty(cpVect, surface_v, SurfaceVelocity, cpTrue);
CP_DefineShapeStructProperty(cpDataPointer, data, UserData, cpFalse);
CP_DefineShapeStructProperty(cpCollisionType, collision_type, CollisionType, cpTrue);
CP_DefineShapeStructProperty(cpGroup, group, Group, cpTrue);
CP_DefineShapeStructProperty(cpLayers, layers, Layers, cpTrue);
*/

/// Extended point query info struct. Returned from calling pointQuery on a shape.
var PointQueryExtendedInfo = function(shape)
{
	/// Shape that was hit, NULL if no collision occurred.
	this.shape = shape;
	/// Depth of the point inside the shape.
	this.d = Infinity;
	/// Direction of minimum norm to the shape's surface.
	this.n = vzero;
};

var NearestPointQueryInfo = function(shape, p, d)
{
	/// The nearest shape, NULL if no shape was within range.
	this.shape = shape;
	/// The closest point on the shape's surface. (in world space coordinates)
	this.p = p;
	/// The distance to the point. The distance is negative if the point is inside the shape.
	this.d = d;
};

var SegmentQueryInfo = function(shape, t, n)
{
	/// The shape that was hit, NULL if no collision occured.
	this.shape = shape;
	/// The normalized distance along the query segment in the range [0, 1].
	this.t = t;
	/// The normal of the surface hit.
	this.n = n;
};

/// Get the hit point for a segment query.
SegmentQueryInfo.prototype.hitPoint = function(start, end)
{
	return vlerp(start, end, this.t);
};

/// Get the hit distance for a segment query.
SegmentQueryInfo.prototype.hitDist = function(start, end)
{
	return vdist(start, end) * this.t;
};

// Circles.

var CircleShape = cp.CircleShape = function(body, radius, offset)
{
	this.c = this.tc = offset;
	this.r = radius;
	
	this.type = 'circle';

	Shape.call(this, body);
};

CircleShape.prototype = Object.create(Shape.prototype);

CircleShape.prototype.cacheData = function(p, rot)
{
	//var c = this.tc = vadd(p, vrotate(this.c, rot));
	var c = this.tc = vrotate(this.c, rot).add(p);
	//this.bb = bbNewForCircle(c, this.r);
	var r = this.r;
	this.bb_l = c.x - r;
	this.bb_b = c.y - r;
	this.bb_r = c.x + r;
	this.bb_t = c.y + r;
};

/// Test if a point lies within a shape.
/*CircleShape.prototype.pointQuery = function(p)
{
	var delta = vsub(p, this.tc);
	var distsq = vlengthsq(delta);
	var r = this.r;
	
	if(distsq < r*r){
		var info = new PointQueryExtendedInfo(this);
		
		var dist = Math.sqrt(distsq);
		info.d = r - dist;
		info.n = vmult(delta, 1/dist);
		return info;
	}
};*/

CircleShape.prototype.nearestPointQuery = function(p)
{
	var deltax = p.x - this.tc.x;
	var deltay = p.y - this.tc.y;
	var d = vlength2(deltax, deltay);
	var r = this.r;
	
	var nearestp = new Vect(this.tc.x + deltax * r/d, this.tc.y + deltay * r/d);
	return new NearestPointQueryInfo(this, nearestp, d - r);
};

var circleSegmentQuery = function(shape, center, r, a, b, info)
{
	// offset the line to be relative to the circle
	a = vsub(a, center);
	b = vsub(b, center);
	
	var qa = vdot(a, a) - 2*vdot(a, b) + vdot(b, b);
	var qb = -2*vdot(a, a) + 2*vdot(a, b);
	var qc = vdot(a, a) - r*r;
	
	var det = qb*qb - 4*qa*qc;
	
	if(det >= 0)
	{
		var t = (-qb - Math.sqrt(det))/(2*qa);
		if(0 <= t && t <= 1){
			return new SegmentQueryInfo(shape, t, vnormalize(vlerp(a, b, t)));
		}
	}
};

CircleShape.prototype.segmentQuery = function(a, b)
{
	return circleSegmentQuery(this, this.tc, this.r, a, b);
};

// The C API has these, and also getters. Its not idiomatic to
// write getters and setters in JS.
/*
CircleShape.prototype.setRadius = function(radius)
{
	this.r = radius;
}

CircleShape.prototype.setOffset = function(offset)
{
	this.c = offset;
}*/

// Segment shape

var SegmentShape = cp.SegmentShape = function(body, a, b, r)
{
	this.a = a;
	this.b = b;
	this.n = vperp(vnormalize(vsub(b, a)));

	this.ta = this.tb = this.tn = null;
	
	this.r = r;
	
	this.a_tangent = vzero;
	this.b_tangent = vzero;
	
	this.type = 'segment';
	Shape.call(this, body);
};

SegmentShape.prototype = Object.create(Shape.prototype);

SegmentShape.prototype.cacheData = function(p, rot)
{
	this.ta = vadd(p, vrotate(this.a, rot));
	this.tb = vadd(p, vrotate(this.b, rot));
	this.tn = vrotate(this.n, rot);
	
	var l,r,b,t;
	
	if(this.ta.x < this.tb.x){
		l = this.ta.x;
		r = this.tb.x;
	} else {
		l = this.tb.x;
		r = this.ta.x;
	}
	
	if(this.ta.y < this.tb.y){
		b = this.ta.y;
		t = this.tb.y;
	} else {
		b = this.tb.y;
		t = this.ta.y;
	}
	
	var rad = this.r;

	this.bb_l = l - rad;
	this.bb_b = b - rad;
	this.bb_r = r + rad;
	this.bb_t = t + rad;
};

SegmentShape.prototype.nearestPointQuery = function(p)
{
	var closest = closestPointOnSegment(p, this.ta, this.tb);
		
	var deltax = p.x - closest.x;
	var deltay = p.y - closest.y;
	var d = vlength2(deltax, deltay);
	var r = this.r;
	
	var nearestp = (d ? vadd(closest, vmult(new Vect(deltax, deltay), r/d)) : closest);
	return new NearestPointQueryInfo(this, nearestp, d - r);
};

SegmentShape.prototype.segmentQuery = function(a, b)
{
	var n = this.tn;
	var d = vdot(vsub(this.ta, a), n);
	var r = this.r;
	
	var flipped_n = (d > 0 ? vneg(n) : n);
	var n_offset = vsub(vmult(flipped_n, r), a);
	
	var seg_a = vadd(this.ta, n_offset);
	var seg_b = vadd(this.tb, n_offset);
	var delta = vsub(b, a);
	
	if(vcross(delta, seg_a)*vcross(delta, seg_b) <= 0){
		var d_offset = d + (d > 0 ? -r : r);
		var ad = -d_offset;
		var bd = vdot(delta, n) - d_offset;
		
		if(ad*bd < 0){
			return new SegmentQueryInfo(this, ad/(ad - bd), flipped_n);
		}
	} else if(r !== 0){
		var info1 = circleSegmentQuery(this, this.ta, this.r, a, b);
		var info2 = circleSegmentQuery(this, this.tb, this.r, a, b);
		
		if (info1){
			return info2 && info2.t < info1.t ? info2 : info1;
		} else {
			return info2;
		}
	}
};

SegmentShape.prototype.setNeighbors = function(prev, next)
{
	this.a_tangent = vsub(prev, this.a);
	this.b_tangent = vsub(next, this.b);
};

SegmentShape.prototype.setEndpoints = function(a, b)
{
	this.a = a;
	this.b = b;
	this.n = vperp(vnormalize(vsub(b, a)));
};

/*
cpSegmentShapeSetRadius(cpShape *shape, cpFloat radius)
{
	this.r = radius;
}*/

/*
CP_DeclareShapeGetter(cpSegmentShape, cpVect, A);
CP_DeclareShapeGetter(cpSegmentShape, cpVect, B);
CP_DeclareShapeGetter(cpSegmentShape, cpVect, Normal);
CP_DeclareShapeGetter(cpSegmentShape, cpFloat, Radius);
*/

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
 
/// Check that a set of vertexes is convex and has a clockwise winding.
var polyValidate = function(verts)
{
	var len = verts.length;
	for(var i=0; i<len; i+=2){
		var ax = verts[i];
	 	var ay = verts[i+1];
		var bx = verts[(i+2)%len];
		var by = verts[(i+3)%len];
		var cx = verts[(i+4)%len];
		var cy = verts[(i+5)%len];
		
		//if(vcross(vsub(b, a), vsub(c, b)) > 0){
		if(vcross2(bx - ax, by - ay, cx - bx, cy - by) > 0){
			return false;
		}
	}
	
	return true;
};

/// Initialize a polygon shape.
/// The vertexes must be convex and have a clockwise winding.
var PolyShape = cp.PolyShape = function(body, verts, offset)
{
	this.setVerts(verts, offset);
	this.type = 'poly';
	Shape.call(this, body);
};

PolyShape.prototype = Object.create(Shape.prototype);

var SplittingPlane = function(n, d)
{
	this.n = n;
	this.d = d;
};

SplittingPlane.prototype.compare = function(v)
{
	return vdot(this.n, v) - this.d;
};

PolyShape.prototype.setVerts = function(verts, offset)
{
	assert(verts.length >= 4, "Polygons require some verts");
	assert(typeof(verts[0]) === 'number',
			'Polygon verticies should be specified in a flattened list (eg [x1,y1,x2,y2,x3,y3,...])');

	// Fail if the user attempts to pass a concave poly, or a bad winding.
	assert(polyValidate(verts), "Polygon is concave or has a reversed winding. Consider using cpConvexHull()");
	
	var len = verts.length;
	var numVerts = len >> 1;

	// This a pretty bad way to do this in javascript. As a first pass, I want to keep
	// the code similar to the C.
	this.verts = new Array(len);
	this.tVerts = new Array(len);
	this.planes = new Array(numVerts);
	this.tPlanes = new Array(numVerts);
	
	for(var i=0; i<len; i+=2){
		//var a = vadd(offset, verts[i]);
		//var b = vadd(offset, verts[(i+1)%numVerts]);
		var ax = verts[i] + offset.x;
	 	var ay = verts[i+1] + offset.y;
		var bx = verts[(i+2)%len] + offset.x;
		var by = verts[(i+3)%len] + offset.y;

		// Inefficient, but only called during object initialization.
		var n = vnormalize(vperp(new Vect(bx-ax, by-ay)));

		this.verts[i  ] = ax;
		this.verts[i+1] = ay;
		this.planes[i>>1] = new SplittingPlane(n, vdot2(n.x, n.y, ax, ay));
		this.tPlanes[i>>1] = new SplittingPlane(new Vect(0,0), 0);
	}
};

/// Initialize a box shaped polygon shape.
var BoxShape = cp.BoxShape = function(body, width, height)
{
	var hw = width/2;
	var hh = height/2;
	
	return BoxShape2(body, new BB(-hw, -hh, hw, hh));
};

/// Initialize an offset box shaped polygon shape.
var BoxShape2 = cp.BoxShape2 = function(body, box)
{
	var verts = [
		box.l, box.b,
		box.l, box.t,
		box.r, box.t,
		box.r, box.b,
	];
	
	return new PolyShape(body, verts, vzero);
};

PolyShape.prototype.transformVerts = function(p, rot)
{
	var src = this.verts;
	var dst = this.tVerts;
	
	var l = Infinity, r = -Infinity;
	var b = Infinity, t = -Infinity;
	
	for(var i=0; i<src.length; i+=2){
		//var v = vadd(p, vrotate(src[i], rot));
		var x = src[i];
	 	var y = src[i+1];

		var vx = p.x + x*rot.x - y*rot.y;
		var vy = p.y + x*rot.y + y*rot.x;

		//console.log('(' + x + ',' + y + ') -> (' + vx + ',' + vy + ')');
		
		dst[i] = vx;
		dst[i+1] = vy;

		l = min(l, vx);
		r = max(r, vx);
		b = min(b, vy);
		t = max(t, vy);
	}

	this.bb_l = l;
	this.bb_b = b;
	this.bb_r = r;
	this.bb_t = t;
};

PolyShape.prototype.transformAxes = function(p, rot)
{
	var src = this.planes;
	var dst = this.tPlanes;
	
	for(var i=0; i<src.length; i++){
		var n = vrotate(src[i].n, rot);
		dst[i].n = n;
		dst[i].d = vdot(p, n) + src[i].d;
	}
};

PolyShape.prototype.cacheData = function(p, rot)
{
	this.transformAxes(p, rot);
	this.transformVerts(p, rot);
};

PolyShape.prototype.nearestPointQuery = function(p)
{
	var planes = this.tPlanes;
	var verts = this.tVerts;
	
	var v0x = verts[verts.length - 2];
	var v0y = verts[verts.length - 1];
	var minDist = Infinity;
	var closestPoint = vzero;
	var outside = false;
	
	for(var i=0; i<planes.length; i++){
		if(planes[i].compare(p) > 0) outside = true;
		
		var v1x = verts[i*2];
		var v1y = verts[i*2 + 1];
		var closest = closestPointOnSegment2(p.x, p.y, v0x, v0y, v1x, v1y);
		
		var dist = vdist(p, closest);
		if(dist < minDist){
			minDist = dist;
			closestPoint = closest;
		}
		
		v0x = v1x;
		v0y = v1y;
	}
	
	return new NearestPointQueryInfo(this, closestPoint, (outside ? minDist : -minDist));
};

PolyShape.prototype.segmentQuery = function(a, b)
{
	var axes = this.tPlanes;
	var verts = this.tVerts;
	var numVerts = axes.length;
	var len = numVerts * 2;
	
	for(var i=0; i<numVerts; i++){
		var n = axes[i].n;
		var an = vdot(a, n);
		if(axes[i].d > an) continue;
		
		var bn = vdot(b, n);
		var t = (axes[i].d - an)/(bn - an);
		if(t < 0 || 1 < t) continue;
		
		var point = vlerp(a, b, t);
		var dt = -vcross(n, point);
		var dtMin = -vcross2(n.x, n.y, verts[i*2], verts[i*2+1]);
		var dtMax = -vcross2(n.x, n.y, verts[(i*2+2)%len], verts[(i*2+3)%len]);

		if(dtMin <= dt && dt <= dtMax){
			// josephg: In the original C code, this function keeps
			// looping through axes after finding a match. I *think*
			// this code is equivalent...
			return new SegmentQueryInfo(this, t, n);
		}
	}
};

PolyShape.prototype.valueOnAxis = function(n, d)
{
	var verts = this.tVerts;
	var m = vdot2(n.x, n.y, verts[0], verts[1]);
	
	for(var i=2; i<verts.length; i+=2){
		m = min(m, vdot2(n.x, n.y, verts[i], verts[i+1]));
	}
	
	return m - d;
};

PolyShape.prototype.containsVert = function(vx, vy)
{
	var planes = this.tPlanes;
	
	for(var i=0; i<planes.length; i++){
		var n = planes[i].n;
		var dist = vdot2(n.x, n.y, vx, vy) - planes[i].d;
		if(dist > 0) return false;
	}
	
	return true;
};

PolyShape.prototype.containsVertPartial = function(vx, vy, n)
{
	var planes = this.tPlanes;
	
	for(var i=0; i<planes.length; i++){
		var n2 = planes[i].n;
		if(vdot(n2, n) < 0) continue;
		var dist = vdot2(n2.x, n2.y, vx, vy) - planes[i].d;
		if(dist > 0) return false;
	}
	
	return true;
};

// These methods are provided for API compatibility with Chipmunk. I recommend against using
// them - just access the poly.verts list directly.
PolyShape.prototype.getNumVerts = function() { return this.verts.length / 2; };
PolyShape.prototype.getVert = function(i)
{
	return new Vect(this.verts[i * 2], this.verts[i * 2 + 1]);
};

/* Copyright (c) 2007 Scott Lembcke
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/// @defgroup cpBody cpBody
/// Chipmunk's rigid body type. Rigid bodies hold the physical properties of an object like
/// it's mass, and position and velocity of it's center of gravity. They don't have an shape on their own.
/// They are given a shape by creating collision shapes (cpShape) that point to the body.
/// @{

var Body = cp.Body = function(m, i) {
	/// Mass of the body.
	/// Must agree with cpBody.m_inv! Use body.setMass() when changing the mass for this reason.
	//this.m;
	/// Mass inverse.
	//this.m_inv;

	/// Moment of inertia of the body.
	/// Must agree with cpBody.i_inv! Use body.setMoment() when changing the moment for this reason.
	//this.i;
	/// Moment of inertia inverse.
	//this.i_inv;

	/// Position of the rigid body's center of gravity.
	this.p = new Vect(0,0);
	/// Velocity of the rigid body's center of gravity.
	this.vx = this.vy = 0;
	/// Force acting on the rigid body's center of gravity.
	this.f = new Vect(0,0);

	/// Rotation of the body around it's center of gravity in radians.
	/// Must agree with cpBody.rot! Use cpBodySetAngle() when changing the angle for this reason.
	//this.a;
	/// Angular velocity of the body around it's center of gravity in radians/second.
	this.w = 0;
	/// Torque applied to the body around it's center of gravity.
	this.t = 0;

	/// Cached unit length vector representing the angle of the body.
	/// Used for fast rotations using cpvrotate().
	//cpVect rot;

	/// Maximum velocity allowed when updating the velocity.
	this.v_limit = Infinity;
	/// Maximum rotational rate (in radians/second) allowed when updating the angular velocity.
	this.w_limit = Infinity;

	// This stuff is all private.
	this.v_biasx = this.v_biasy = 0;
	this.w_bias = 0;

	this.space = null;

	this.shapeList = [];
	this.arbiterList = null; // These are both wacky linked lists.
	this.constraintList = null;

	// This stuff is used to track information on the collision graph.
	this.nodeRoot = null;
	this.nodeNext = null;
	this.nodeIdleTime = 0;

	// Set this.m and this.m_inv
	this.setMass(m);

	// Set this.i and this.i_inv
	this.setMoment(i);

	// Set this.a and this.rot
	this.rot = new Vect(0,0);
	this.setAngle(0);
};

// I wonder if this should use the constructor style like Body...
var createStaticBody = function()
{
	body = new Body(Infinity, Infinity);
	body.nodeIdleTime = Infinity;

	return body;
};

if (typeof DEBUG !== 'undefined' && DEBUG) {
	var v_assert_nan = function(v, message){assert(v.x == v.x && v.y == v.y, message); };
	var v_assert_infinite = function(v, message){assert(Math.abs(v.x) !== Infinity && Math.abs(v.y) !== Infinity, message);};
	var v_assert_sane = function(v, message){v_assert_nan(v, message); v_assert_infinite(v, message);};

	Body.prototype.sanityCheck = function()
	{
		assert(this.m === this.m && this.m_inv === this.m_inv, "Body's mass is invalid.");
		assert(this.i === this.i && this.i_inv === this.i_inv, "Body's moment is invalid.");

		v_assert_sane(this.p, "Body's position is invalid.");
		v_assert_sane(this.f, "Body's force is invalid.");
		assert(this.vx === this.vx && Math.abs(this.vx) !== Infinity, "Body's velocity is invalid.");
		assert(this.vy === this.vy && Math.abs(this.vy) !== Infinity, "Body's velocity is invalid.");

		assert(this.a === this.a && Math.abs(this.a) !== Infinity, "Body's angle is invalid.");
		assert(this.w === this.w && Math.abs(this.w) !== Infinity, "Body's angular velocity is invalid.");
		assert(this.t === this.t && Math.abs(this.t) !== Infinity, "Body's torque is invalid.");

		v_assert_sane(this.rot, "Body's rotation vector is invalid.");

		assert(this.v_limit === this.v_limit, "Body's velocity limit is invalid.");
		assert(this.w_limit === this.w_limit, "Body's angular velocity limit is invalid.");
	};
} else {
	Body.prototype.sanityCheck = function(){};
}

Body.prototype.getPos = function() { return this.p; };
Body.prototype.getVel = function() { return new Vect(this.vx, this.vy); };
Body.prototype.getAngVel = function() { return this.w; };

/// Returns true if the body is sleeping.
Body.prototype.isSleeping = function()
{
	return this.nodeRoot !== null;
};

/// Returns true if the body is static.
Body.prototype.isStatic = function()
{
	return this.nodeIdleTime === Infinity;
};

/// Returns true if the body has not been added to a space.
Body.prototype.isRogue = function()
{
	return this.space === null;
};

// It would be nicer to use defineProperty for this, but its about 30x slower:
// http://jsperf.com/defineproperty-vs-setter
Body.prototype.setMass = function(mass)
{
	assert(mass > 0, "Mass must be positive and non-zero.");

	//activate is defined in cpSpaceComponent
	this.activate();
	this.m = mass;
	this.m_inv = 1/mass;
};

Body.prototype.setMoment = function(moment)
{
	assert(moment > 0, "Moment of Inertia must be positive and non-zero.");

	this.activate();
	this.i = moment;
	this.i_inv = 1/moment;
};

Body.prototype.addShape = function(shape)
{
	this.shapeList.push(shape);
};

Body.prototype.removeShape = function(shape)
{
	// This implementation has a linear time complexity with the number of shapes.
	// The original implementation used linked lists instead, which might be faster if
	// you're constantly editing the shape of a body. I expect most bodies will never
	// have their shape edited, so I'm just going to use the simplest possible implemention.
	deleteObjFromList(this.shapeList, shape);
};

var filterConstraints = function(node, body, filter)
{
	if(node === filter){
		return node.next(body);
	} else if(node.a === body){
		node.next_a = filterConstraints(node.next_a, body, filter);
	} else {
		node.next_b = filterConstraints(node.next_b, body, filter);
	}

	return node;
};

Body.prototype.removeConstraint = function(constraint)
{
	// The constraint must be in the constraints list when this is called.
	this.constraintList = filterConstraints(this.constraintList, this, constraint);
};

Body.prototype.setPos = function(pos)
{
	this.activate();
	this.sanityCheck();
	// If I allow the position to be set to vzero, vzero will get changed.
	if (pos === vzero) {
		pos = cp.v(0,0);
	}
	this.p = pos;
};

Body.prototype.setVel = function(velocity)
{
	this.activate();
	this.vx = velocity.x;
	this.vy = velocity.y;
};

Body.prototype.setAngVel = function(w)
{
	this.activate();
	this.w = w;
};

Body.prototype.setAngleInternal = function(angle)
{
	assert(!isNaN(angle), "Internal Error: Attempting to set body's angle to NaN");
	this.a = angle;//fmod(a, (cpFloat)M_PI*2.0f);

	//this.rot = vforangle(angle);
	this.rot.x = Math.cos(angle);
	this.rot.y = Math.sin(angle);
};

Body.prototype.setAngle = function(angle)
{
	this.activate();
	this.sanityCheck();
	this.setAngleInternal(angle);
};

Body.prototype.velocity_func = function(gravity, damping, dt)
{
	//this.v = vclamp(vadd(vmult(this.v, damping), vmult(vadd(gravity, vmult(this.f, this.m_inv)), dt)), this.v_limit);
	var vx = this.vx * damping + (gravity.x + this.f.x * this.m_inv) * dt;
	var vy = this.vy * damping + (gravity.y + this.f.y * this.m_inv) * dt;

	//var v = vclamp(new Vect(vx, vy), this.v_limit);
	//this.vx = v.x; this.vy = v.y;
	var v_limit = this.v_limit;
	var lensq = vx * vx + vy * vy;
	var scale = (lensq > v_limit*v_limit) ? v_limit / Math.sqrt(lensq) : 1;
	this.vx = vx * scale;
	this.vy = vy * scale;

	var w_limit = this.w_limit;
	this.w = clamp(this.w*damping + this.t*this.i_inv*dt, -w_limit, w_limit);

	this.sanityCheck();
};

Body.prototype.position_func = function(dt)
{
	//this.p = vadd(this.p, vmult(vadd(this.v, this.v_bias), dt));

	//this.p = this.p + (this.v + this.v_bias) * dt;
	this.p.x += (this.vx + this.v_biasx) * dt;
	this.p.y += (this.vy + this.v_biasy) * dt;

	this.setAngleInternal(this.a + (this.w + this.w_bias)*dt);

	this.v_biasx = this.v_biasy = 0;
	this.w_bias = 0;

	this.sanityCheck();
};

Body.prototype.resetForces = function()
{
	this.activate();
	this.f = new Vect(0,0);
	this.t = 0;
};

Body.prototype.applyForce = function(force, r)
{
	this.activate();
	this.f = vadd(this.f, force);
	this.t += vcross(r, force);
};

Body.prototype.applyImpulse = function(j, r)
{
	this.activate();
	apply_impulse(this, j.x, j.y, r);
};

Body.prototype.getVelAtPoint = function(r)
{
	return vadd(new Vect(this.vx, this.vy), vmult(vperp(r), this.w));
};

/// Get the velocity on a body (in world units) at a point on the body in world coordinates.
Body.prototype.getVelAtWorldPoint = function(point)
{
	return this.getVelAtPoint(vsub(point, this.p));
};

/// Get the velocity on a body (in world units) at a point on the body in local coordinates.
Body.prototype.getVelAtLocalPoint = function(point)
{
	return this.getVelAtPoint(vrotate(point, this.rot));
};

Body.prototype.eachShape = function(func)
{
	for(var i = 0, len = this.shapeList.length; i < len; i++) {
		func(this.shapeList[i]);
	}
};

Body.prototype.eachConstraint = function(func)
{
	var constraint = this.constraintList;
	while(constraint) {
		var next = constraint.next(this);
		func(constraint);
		constraint = next;
	}
};

Body.prototype.eachArbiter = function(func)
{
	var arb = this.arbiterList;
	while(arb){
		var next = arb.next(this);

		arb.swappedColl = (this === arb.body_b);
		func(arb);

		arb = next;
	}
};

/// Convert body relative/local coordinates to absolute/world coordinates.
Body.prototype.local2World = function(v)
{
	return vadd(this.p, vrotate(v, this.rot));
};

/// Convert body absolute/world coordinates to	relative/local coordinates.
Body.prototype.world2Local = function(v)
{
	return vunrotate(vsub(v, this.p), this.rot);
};

/// Get the kinetic energy of a body.
Body.prototype.kineticEnergy = function()
{
	// Need to do some fudging to avoid NaNs
	var vsq = this.vx*this.vx + this.vy*this.vy;
	var wsq = this.w * this.w;
	return (vsq ? vsq*this.m : 0) + (wsq ? wsq*this.i : 0);
};

/* Copyright (c) 2010 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
	@defgroup cpSpatialIndex cpSpatialIndex
	
	Spatial indexes are data structures that are used to accelerate collision detection
	and spatial queries. Chipmunk provides a number of spatial index algorithms to pick from
	and they are programmed in a generic way so that you can use them for holding more than
	just Shapes.
	
	It works by using pointers to the objects you add and using a callback to ask your code
	for bounding boxes when it needs them. Several types of queries can be performed an index as well
	as reindexing and full collision information. All communication to the spatial indexes is performed
	through callback functions.
	
	Spatial indexes should be treated as opaque structs.
	This means you shouldn't be reading any of the fields directly.

	All spatial indexes define the following methods:
		
	// The number of objects in the spatial index.
	count = 0;

	// Iterate the objects in the spatial index. @c func will be called once for each object.
	each(func);
	
	// Returns true if the spatial index contains the given object.
	// Most spatial indexes use hashed storage, so you must provide a hash value too.
	contains(obj, hashid);

	// Add an object to a spatial index.
	insert(obj, hashid);

	// Remove an object from a spatial index.
	remove(obj, hashid);
	
	// Perform a full reindex of a spatial index.
	reindex();

	// Reindex a single object in the spatial index.
	reindexObject(obj, hashid);

	// Perform a point query against the spatial index, calling @c func for each potential match.
	// A pointer to the point will be passed as @c obj1 of @c func.
	// func(shape);
	pointQuery(point, func);

	// Perform a segment query against the spatial index, calling @c func for each potential match.
	// func(shape);
	segmentQuery(vect a, vect b, t_exit, func);

	// Perform a rectangle query against the spatial index, calling @c func for each potential match.
	// func(shape);
	query(bb, func);

	// Simultaneously reindex and find all colliding objects.
	// @c func will be called once for each potentially overlapping pair of objects found.
	// If the spatial index was initialized with a static index, it will collide it's objects against that as well.
	reindexQuery(func);
*/

var SpatialIndex = cp.SpatialIndex = function(staticIndex)
{
	this.staticIndex = staticIndex;
	

	if(staticIndex){
		if(staticIndex.dynamicIndex){
			throw new Error("This static index is already associated with a dynamic index.");
		}
		staticIndex.dynamicIndex = this;
	}
};

// Collide the objects in an index against the objects in a staticIndex using the query callback function.
SpatialIndex.prototype.collideStatic = function(staticIndex, func)
{
	if(staticIndex.count > 0){
		var query = staticIndex.query;

		this.each(function(obj) {
			query(obj, new BB(obj.bb_l, obj.bb_b, obj.bb_r, obj.bb_t), func);
		});
	}
};


/* Copyright (c) 2009 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// This file implements a modified AABB tree for collision detection.

var BBTree = cp.BBTree = function(staticIndex)
{
	SpatialIndex.call(this, staticIndex);
	
	this.velocityFunc = null;

	// This is a hash from object ID -> object for the objects stored in the BBTree.
	this.leaves = {};
	// A count of the number of leaves in the BBTree.
	this.count = 0;

	this.root = null;
	
	// A linked list containing an object pool of tree nodes and pairs.
	this.pooledNodes = null;
	this.pooledPairs = null;
	
	this.stamp = 0;
};

BBTree.prototype = Object.create(SpatialIndex.prototype);

var numNodes = 0;

var Node = function(tree, a, b)
{
	this.obj = null;
	this.bb_l = min(a.bb_l, b.bb_l);
	this.bb_b = min(a.bb_b, b.bb_b);
	this.bb_r = max(a.bb_r, b.bb_r);
	this.bb_t = max(a.bb_t, b.bb_t);
	this.parent = null;
	
	this.setA(a);
	this.setB(b);
};

BBTree.prototype.makeNode = function(a, b)
{
	var node = this.pooledNodes;
	if(node){
		this.pooledNodes = node.parent;
		node.constructor(this, a, b);
		return node;
	} else {
		numNodes++;
		return new Node(this, a, b);
	}
};

var numLeaves = 0;
var Leaf = function(tree, obj)
{
	this.obj = obj;
	tree.getBB(obj, this);

	this.parent = null;

	this.stamp = 1;
	this.pairs = null;
	numLeaves++;
};

// **** Misc Functions

BBTree.prototype.getBB = function(obj, dest)
{
	var velocityFunc = this.velocityFunc;
	if(velocityFunc){
		var coef = 0.1;
		var x = (obj.bb_r - obj.bb_l)*coef;
		var y = (obj.bb_t - obj.bb_b)*coef;
		
		var v = vmult(velocityFunc(obj), 0.1);

		dest.bb_l = obj.bb_l + min(-x, v.x);
		dest.bb_b = obj.bb_b + min(-y, v.y);
		dest.bb_r = obj.bb_r + max( x, v.x);
		dest.bb_t = obj.bb_t + max( y, v.y);
	} else {
		dest.bb_l = obj.bb_l;
		dest.bb_b = obj.bb_b;
		dest.bb_r = obj.bb_r;
		dest.bb_t = obj.bb_t;
	}
};

BBTree.prototype.getStamp = function()
{
	var dynamic = this.dynamicIndex;
	return (dynamic && dynamic.stamp ? dynamic.stamp : this.stamp);
};

BBTree.prototype.incrementStamp = function()
{
	if(this.dynamicIndex && this.dynamicIndex.stamp){
		this.dynamicIndex.stamp++;
	} else {
		this.stamp++;
	}
}

// **** Pair/Thread Functions

var numPairs = 0;
// Objects created with constructors are faster than object literals. :(
var Pair = function(leafA, nextA, leafB, nextB)
{
	this.prevA = null;
	this.leafA = leafA;
	this.nextA = nextA;

	this.prevB = null;
	this.leafB = leafB;
	this.nextB = nextB;
};

BBTree.prototype.makePair = function(leafA, nextA, leafB, nextB)
{
	//return new Pair(leafA, nextA, leafB, nextB);
	var pair = this.pooledPairs;
	if (pair)
	{
		this.pooledPairs = pair.prevA;

		pair.prevA = null;
		pair.leafA = leafA;
		pair.nextA = nextA;

		pair.prevB = null;
		pair.leafB = leafB;
		pair.nextB = nextB;

		//pair.constructor(leafA, nextA, leafB, nextB);
		return pair;
	} else {
		numPairs++;
		return new Pair(leafA, nextA, leafB, nextB);
	}
};

Pair.prototype.recycle = function(tree)
{
	this.prevA = tree.pooledPairs;
	tree.pooledPairs = this;
};

var unlinkThread = function(prev, leaf, next)
{
	if(next){
		if(next.leafA === leaf) next.prevA = prev; else next.prevB = prev;
	}
	
	if(prev){
		if(prev.leafA === leaf) prev.nextA = next; else prev.nextB = next;
	} else {
		leaf.pairs = next;
	}
};

Leaf.prototype.clearPairs = function(tree)
{
	var pair = this.pairs,
		next;

	this.pairs = null;
	
	while(pair){
		if(pair.leafA === this){
			next = pair.nextA;
			unlinkThread(pair.prevB, pair.leafB, pair.nextB);
		} else {
			next = pair.nextB;
			unlinkThread(pair.prevA, pair.leafA, pair.nextA);
		}
		pair.recycle(tree);
		pair = next;
	}
};

var pairInsert = function(a, b, tree)
{
	var nextA = a.pairs, nextB = b.pairs;
	var pair = tree.makePair(a, nextA, b, nextB);
	a.pairs = b.pairs = pair;

	if(nextA){
		if(nextA.leafA === a) nextA.prevA = pair; else nextA.prevB = pair;
	}
	
	if(nextB){
		if(nextB.leafA === b) nextB.prevA = pair; else nextB.prevB = pair;
	}
};

// **** Node Functions

Node.prototype.recycle = function(tree)
{
	this.parent = tree.pooledNodes;
	tree.pooledNodes = this;
};

Leaf.prototype.recycle = function(tree)
{
	// Its not worth the overhead to recycle leaves.
};

Node.prototype.setA = function(value)
{
	this.A = value;
	value.parent = this;
};

Node.prototype.setB = function(value)
{
	this.B = value;
	value.parent = this;
};

Leaf.prototype.isLeaf = true;
Node.prototype.isLeaf = false;

Node.prototype.otherChild = function(child)
{
	return (this.A == child ? this.B : this.A);
};

Node.prototype.replaceChild = function(child, value, tree)
{
	assertSoft(child == this.A || child == this.B, "Node is not a child of parent.");
	
	if(this.A == child){
		this.A.recycle(tree);
		this.setA(value);
	} else {
		this.B.recycle(tree);
		this.setB(value);
	}
	
	for(var node=this; node; node = node.parent){
		//node.bb = bbMerge(node.A.bb, node.B.bb);
		var a = node.A;
		var b = node.B;
		node.bb_l = min(a.bb_l, b.bb_l);
		node.bb_b = min(a.bb_b, b.bb_b);
		node.bb_r = max(a.bb_r, b.bb_r);
		node.bb_t = max(a.bb_t, b.bb_t);
	}
};

Node.prototype.bbArea = Leaf.prototype.bbArea = function()
{
	return (this.bb_r - this.bb_l)*(this.bb_t - this.bb_b);
};

var bbTreeMergedArea = function(a, b)
{
	return (max(a.bb_r, b.bb_r) - min(a.bb_l, b.bb_l))*(max(a.bb_t, b.bb_t) - min(a.bb_b, b.bb_b));
};

// **** Subtree Functions

// Would it be better to make these functions instance methods on Node and Leaf?

var bbProximity = function(a, b)
{
	return Math.abs(a.bb_l + a.bb_r - b.bb_l - b.bb_r) + Math.abs(a.bb_b + a.bb_t - b.bb_b - b.bb_t);
};

var subtreeInsert = function(subtree, leaf, tree)
{
//	var s = new Error().stack;
//	traces[s] = traces[s] ? traces[s]+1 : 1;

	if(subtree == null){
		return leaf;
	} else if(subtree.isLeaf){
		return tree.makeNode(leaf, subtree);
	} else {
		var cost_a = subtree.B.bbArea() + bbTreeMergedArea(subtree.A, leaf);
		var cost_b = subtree.A.bbArea() + bbTreeMergedArea(subtree.B, leaf);
		
		if(cost_a === cost_b){
			cost_a = bbProximity(subtree.A, leaf);
			cost_b = bbProximity(subtree.B, leaf);
		}	

		if(cost_b < cost_a){
			subtree.setB(subtreeInsert(subtree.B, leaf, tree));
		} else {
			subtree.setA(subtreeInsert(subtree.A, leaf, tree));
		}
		
//		subtree.bb = bbMerge(subtree.bb, leaf.bb);
		subtree.bb_l = min(subtree.bb_l, leaf.bb_l);
		subtree.bb_b = min(subtree.bb_b, leaf.bb_b);
		subtree.bb_r = max(subtree.bb_r, leaf.bb_r);
		subtree.bb_t = max(subtree.bb_t, leaf.bb_t);

		return subtree;
	}
};

Node.prototype.intersectsBB = Leaf.prototype.intersectsBB = function(bb)
{
	return (this.bb_l <= bb.r && bb.l <= this.bb_r && this.bb_b <= bb.t && bb.b <= this.bb_t);
};

var subtreeQuery = function(subtree, bb, func)
{
	//if(bbIntersectsBB(subtree.bb, bb)){
	if(subtree.intersectsBB(bb)){
		if(subtree.isLeaf){
			func(subtree.obj);
		} else {
			subtreeQuery(subtree.A, bb, func);
			subtreeQuery(subtree.B, bb, func);
		}
	}
};

/// Returns the fraction along the segment query the node hits. Returns Infinity if it doesn't hit.
var nodeSegmentQuery = function(node, a, b)
{
	var idx = 1/(b.x - a.x);
	var tx1 = (node.bb_l == a.x ? -Infinity : (node.bb_l - a.x)*idx);
	var tx2 = (node.bb_r == a.x ?  Infinity : (node.bb_r - a.x)*idx);
	var txmin = min(tx1, tx2);
	var txmax = max(tx1, tx2);
	
	var idy = 1/(b.y - a.y);
	var ty1 = (node.bb_b == a.y ? -Infinity : (node.bb_b - a.y)*idy);
	var ty2 = (node.bb_t == a.y ?  Infinity : (node.bb_t - a.y)*idy);
	var tymin = min(ty1, ty2);
	var tymax = max(ty1, ty2);
	
	if(tymin <= txmax && txmin <= tymax){
		var min_ = max(txmin, tymin);
		var max_ = min(txmax, tymax);
		
		if(0.0 <= max_ && min_ <= 1.0) return max(min_, 0.0);
	}
	
	return Infinity;
};

var subtreeSegmentQuery = function(subtree, a, b, t_exit, func)
{
	if(subtree.isLeaf){
		return func(subtree.obj);
	} else {
		var t_a = nodeSegmentQuery(subtree.A, a, b);
		var t_b = nodeSegmentQuery(subtree.B, a, b);
		
		if(t_a < t_b){
			if(t_a < t_exit) t_exit = min(t_exit, subtreeSegmentQuery(subtree.A, a, b, t_exit, func));
			if(t_b < t_exit) t_exit = min(t_exit, subtreeSegmentQuery(subtree.B, a, b, t_exit, func));
		} else {
			if(t_b < t_exit) t_exit = min(t_exit, subtreeSegmentQuery(subtree.B, a, b, t_exit, func));
			if(t_a < t_exit) t_exit = min(t_exit, subtreeSegmentQuery(subtree.A, a, b, t_exit, func));
		}
		
		return t_exit;
	}
};

BBTree.prototype.subtreeRecycle = function(node)
{
	if(node.isLeaf){
		this.subtreeRecycle(node.A);
		this.subtreeRecycle(node.B);
		node.recycle(this);
	}
};

var subtreeRemove = function(subtree, leaf, tree)
{
	if(leaf == subtree){
		return null;
	} else {
		var parent = leaf.parent;
		if(parent == subtree){
			var other = subtree.otherChild(leaf);
			other.parent = subtree.parent;
			subtree.recycle(tree);
			return other;
		} else {
			parent.parent.replaceChild(parent, parent.otherChild(leaf), tree);
			return subtree;
		}
	}
};

// **** Marking Functions

/*
typedef struct MarkContext {
	bbTree *tree;
	Node *staticRoot;
	cpSpatialIndexQueryFunc func;
} MarkContext;
*/

var bbTreeIntersectsNode = function(a, b)
{
	return (a.bb_l <= b.bb_r && b.bb_l <= a.bb_r && a.bb_b <= b.bb_t && b.bb_b <= a.bb_t);
};

Leaf.prototype.markLeafQuery = function(leaf, left, tree, func)
{
	if(bbTreeIntersectsNode(leaf, this)){
    if(left){
      pairInsert(leaf, this, tree);
    } else {
      if(this.stamp < leaf.stamp) pairInsert(this, leaf, tree);
      if(func) func(leaf.obj, this.obj);
    }
  }
};

Node.prototype.markLeafQuery = function(leaf, left, tree, func)
{
	if(bbTreeIntersectsNode(leaf, this)){
    this.A.markLeafQuery(leaf, left, tree, func);
    this.B.markLeafQuery(leaf, left, tree, func);
	}
};

Leaf.prototype.markSubtree = function(tree, staticRoot, func)
{
	if(this.stamp == tree.getStamp()){
		if(staticRoot) staticRoot.markLeafQuery(this, false, tree, func);
		
		for(var node = this; node.parent; node = node.parent){
			if(node == node.parent.A){
				node.parent.B.markLeafQuery(this, true, tree, func);
			} else {
				node.parent.A.markLeafQuery(this, false, tree, func);
			}
		}
	} else {
		var pair = this.pairs;
		while(pair){
			if(this === pair.leafB){
				if(func) func(pair.leafA.obj, this.obj);
				pair = pair.nextB;
			} else {
				pair = pair.nextA;
			}
		}
	}
};

Node.prototype.markSubtree = function(tree, staticRoot, func)
{
  this.A.markSubtree(tree, staticRoot, func);
  this.B.markSubtree(tree, staticRoot, func);
};

// **** Leaf Functions

Leaf.prototype.containsObj = function(obj)
{
	return (this.bb_l <= obj.bb_l && this.bb_r >= obj.bb_r && this.bb_b <= obj.bb_b && this.bb_t >= obj.bb_t);
};

Leaf.prototype.update = function(tree)
{
	var root = tree.root;
	var obj = this.obj;

	//if(!bbContainsBB(this.bb, bb)){
	if(!this.containsObj(obj)){
		tree.getBB(this.obj, this);
		
		root = subtreeRemove(root, this, tree);
		tree.root = subtreeInsert(root, this, tree);
		
		this.clearPairs(tree);
		this.stamp = tree.getStamp();
		
		return true;
	}
	
	return false;
};

Leaf.prototype.addPairs = function(tree)
{
	var dynamicIndex = tree.dynamicIndex;
	if(dynamicIndex){
		var dynamicRoot = dynamicIndex.root;
		if(dynamicRoot){
			dynamicRoot.markLeafQuery(this, true, dynamicIndex, null);
		}
	} else {
		var staticRoot = tree.staticIndex.root;
		this.markSubtree(tree, staticRoot, null);
	}
};

// **** Insert/Remove

BBTree.prototype.insert = function(obj, hashid)
{
	var leaf = new Leaf(this, obj);

	this.leaves[hashid] = leaf;
	this.root = subtreeInsert(this.root, leaf, this);
	this.count++;
	
	leaf.stamp = this.getStamp();
	leaf.addPairs(this);
	this.incrementStamp();
};

BBTree.prototype.remove = function(obj, hashid)
{
	var leaf = this.leaves[hashid];

	delete this.leaves[hashid];
	this.root = subtreeRemove(this.root, leaf, this);
	this.count--;

	leaf.clearPairs(this);
	leaf.recycle(this);
};

BBTree.prototype.contains = function(obj, hashid)
{
	return this.leaves[hashid] != null;
};

// **** Reindex
var voidQueryFunc = function(obj1, obj2){};

BBTree.prototype.reindexQuery = function(func)
{
	if(!this.root) return;
	
	// LeafUpdate() may modify this.root. Don't cache it.
	var hashid,
		leaves = this.leaves;
	for (hashid in leaves)
	{
		leaves[hashid].update(this);
	}
	
	var staticIndex = this.staticIndex;
	var staticRoot = staticIndex && staticIndex.root;
	
	this.root.markSubtree(this, staticRoot, func);
	if(staticIndex && !staticRoot) this.collideStatic(this, staticIndex, func);
	
	this.incrementStamp();
};

BBTree.prototype.reindex = function()
{
	this.reindexQuery(voidQueryFunc);
};

BBTree.prototype.reindexObject = function(obj, hashid)
{
	var leaf = this.leaves[hashid];
	if(leaf){
		if(leaf.update(this)) leaf.addPairs(this);
		this.incrementStamp();
	}
};

// **** Query

// This has since been removed from upstream Chipmunk - which recommends you just use query() below
// directly.
BBTree.prototype.pointQuery = function(point, func)
{
	this.query(new BB(point.x, point.y, point.x, point.y), func);
};

BBTree.prototype.segmentQuery = function(a, b, t_exit, func)
{
	if(this.root) subtreeSegmentQuery(this.root, a, b, t_exit, func);
};

BBTree.prototype.query = function(bb, func)
{
	if(this.root) subtreeQuery(this.root, bb, func);
};

// **** Misc

BBTree.prototype.count = function()
{
	return this.count;
};

BBTree.prototype.each = function(func)
{
	var hashid;
	for(hashid in this.leaves)
	{
		func(this.leaves[hashid].obj);
	}
};

// **** Tree Optimization

var bbTreeMergedArea2 = function(node, l, b, r, t)
{
	return (max(node.bb_r, r) - min(node.bb_l, l))*(max(node.bb_t, t) - min(node.bb_b, b));
};

var partitionNodes = function(tree, nodes, offset, count)
{
	if(count == 1){
		return nodes[offset];
	} else if(count == 2) {
		return tree.makeNode(nodes[offset], nodes[offset + 1]);
	}
	
	// Find the AABB for these nodes
	//var bb = nodes[offset].bb;
	var node = nodes[offset];
	var bb_l = node.bb_l,
		bb_b = node.bb_b,
		bb_r = node.bb_r,
		bb_t = node.bb_t;

	var end = offset + count;
	for(var i=offset + 1; i<end; i++){
		//bb = bbMerge(bb, nodes[i].bb);
		node = nodes[i];
		bb_l = min(bb_l, node.bb_l);
		bb_b = min(bb_b, node.bb_b);
		bb_r = max(bb_r, node.bb_r);
		bb_t = max(bb_t, node.bb_t);
	}
	
	// Split it on it's longest axis
	var splitWidth = (bb_r - bb_l > bb_t - bb_b);
	
	// Sort the bounds and use the median as the splitting point
	var bounds = new Array(count*2);
	if(splitWidth){
		for(var i=offset; i<end; i++){
			bounds[2*i + 0] = nodes[i].bb_l;
			bounds[2*i + 1] = nodes[i].bb_r;
		}
	} else {
		for(var i=offset; i<end; i++){
			bounds[2*i + 0] = nodes[i].bb_b;
			bounds[2*i + 1] = nodes[i].bb_t;
		}
	}
	
	bounds.sort(function(a, b) {
		// This might run faster if the function was moved out into the global scope.
		return a - b;
	});
	var split = (bounds[count - 1] + bounds[count])*0.5; // use the median as the split

	// Generate the child BBs
	//var a = bb, b = bb;
	var a_l = bb_l, a_b = bb_b, a_r = bb_r, a_t = bb_t;
	var b_l = bb_l, b_b = bb_b, b_r = bb_r, b_t = bb_t;

	if(splitWidth) a_r = b_l = split; else a_t = b_b = split;
	
	// Partition the nodes
	var right = end;
	for(var left=offset; left < right;){
		var node = nodes[left];
//	if(bbMergedArea(node.bb, b) < bbMergedArea(node.bb, a)){
		if(bbTreeMergedArea2(node, b_l, b_b, b_r, b_t) < bbTreeMergedArea2(node, a_l, a_b, a_r, a_t)){
			right--;
			nodes[left] = nodes[right];
			nodes[right] = node;
		} else {
			left++;
		}
	}
	
	if(right == count){
		var node = null;
		for(var i=offset; i<end; i++) node = subtreeInsert(node, nodes[i], tree);
		return node;
	}
	
	// Recurse and build the node!
	return NodeNew(tree,
		partitionNodes(tree, nodes, offset, right - offset),
		partitionNodes(tree, nodes, right, end - right)
	);
};

//static void
//bbTreeOptimizeIncremental(bbTree *tree, int passes)
//{
//	for(int i=0; i<passes; i++){
//		Node *root = tree.root;
//		Node *node = root;
//		int bit = 0;
//		unsigned int path = tree.opath;
//		
//		while(!NodeIsLeaf(node)){
//			node = (path&(1<<bit) ? node.a : node.b);
//			bit = (bit + 1)&(sizeof(unsigned int)*8 - 1);
//		}
//		
//		root = subtreeRemove(root, node, tree);
//		tree.root = subtreeInsert(root, node, tree);
//	}
//}

BBTree.prototype.optimize = function()
{
	var nodes = new Array(this.count);
	var i = 0;

	for (var hashid in this.leaves)
	{
		nodes[i++] = this.nodes[hashid];
	}
	
	tree.subtreeRecycle(root);
	this.root = partitionNodes(tree, nodes, nodes.length);
};

// **** Debug Draw

var nodeRender = function(node, depth)
{
	if(!node.isLeaf && depth <= 10){
		nodeRender(node.A, depth + 1);
		nodeRender(node.B, depth + 1);
	}
	
	var str = '';
	for(var i = 0; i < depth; i++) {
		str += ' ';
	}

	console.log(str + node.bb_b + ' ' + node.bb_t);
};

BBTree.prototype.log = function(){
	if(this.root) nodeRender(this.root, 0);
};

/*
static void
NodeRender(Node *node, int depth)
{
	if(!NodeIsLeaf(node) && depth <= 10){
		NodeRender(node.a, depth + 1);
		NodeRender(node.b, depth + 1);
	}
	
	bb bb = node.bb;
	
//	GLfloat v = depth/2.0f;	
//	glColor3f(1.0f - v, v, 0.0f);
	glLineWidth(max(5.0f - depth, 1.0f));
	glBegin(GL_LINES); {
		glVertex2f(bb.l, bb.b);
		glVertex2f(bb.l, bb.t);
		
		glVertex2f(bb.l, bb.t);
		glVertex2f(bb.r, bb.t);
		
		glVertex2f(bb.r, bb.t);
		glVertex2f(bb.r, bb.b);
		
		glVertex2f(bb.r, bb.b);
		glVertex2f(bb.l, bb.b);
	}; glEnd();
}

void
bbTreeRenderDebug(cpSpatialIndex *index){
	if(index.klass != &klass){
		cpAssertWarn(false, "Ignoring bbTreeRenderDebug() call to non-tree spatial index.");
		return;
	}
	
	bbTree *tree = (bbTree *)index;
	if(tree.root) NodeRender(tree.root, 0);
}
*/
/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/// @defgroup cpArbiter cpArbiter
/// The cpArbiter struct controls pairs of colliding shapes.
/// They are also used in conjuction with collision handler callbacks
/// allowing you to retrieve information on the collision and control it.


// **** Collision Handlers
//
// Collision handlers are user-defined objects to describe the behaviour of colliding
// objects.
var CollisionHandler = cp.CollisionHandler = function()
{
	// The collision type
	this.a = this.b = 0;
};

/// Collision begin event callback
/// Returning false from a begin callback causes the collision to be ignored until
/// the the separate callback is called when the objects stop colliding.
CollisionHandler.prototype.begin = function(arb, space){return true;};

/// Collision pre-solve event callback
/// Returning false from a pre-step callback causes the collision to be ignored until the next step.
CollisionHandler.prototype.preSolve = function(arb, space){return true;};

/// Collision post-solve event function callback type.
CollisionHandler.prototype.postSolve = function(arb, space){};
/// Collision separate event function callback type.
CollisionHandler.prototype.separate = function(arb, space){};


var CP_MAX_CONTACTS_PER_ARBITER = 4;

// Arbiter states
//
// Arbiter is active and its the first collision.
//	'first coll'
// Arbiter is active and its not the first collision.
//	'normal',
// Collision has been explicitly ignored.
// Either by returning false from a begin collision handler or calling cpArbiterIgnore().
//	'ignore',
// Collison is no longer active. A space will cache an arbiter for up to cpSpace.collisionPersistence more steps.
//	'cached'

/// A colliding pair of shapes.
var Arbiter = function(a, b) {
	/// Calculated value to use for the elasticity coefficient.
	/// Override in a pre-solve collision handler for custom behavior.
	this.e = 0;
	/// Calculated value to use for the friction coefficient.
	/// Override in a pre-solve collision handler for custom behavior.
	this.u = 0;
	/// Calculated value to use for applying surface velocities.
	/// Override in a pre-solve collision handler for custom behavior.
	this.surface_vr = vzero;
	
	this.a = a; this.body_a = a.body;
	this.b = b; this.body_b = b.body;
	
	this.thread_a_next = this.thread_a_prev = null;
	this.thread_b_next = this.thread_b_prev = null;
	
	this.contacts = null;
	
	this.stamp = 0;
	this.handler = null;
	this.swappedColl = false;
	this.state = 'first coll';
};

Arbiter.prototype.getShapes = function()
{
	if (this.swappedColl){
		return [this.b, this.a];
	}else{
		return [this.a, this.b];
	}
}

/// Calculate the total impulse that was applied by this arbiter.
/// This function should only be called from a post-solve, post-step or cpBodyEachArbiter callback.
Arbiter.prototype.totalImpulse = function()
{
	var contacts = this.contacts;
	var sum = new Vect(0,0);
	
	for(var i=0, count=contacts.length; i<count; i++){
		var con = contacts[i];
		sum.add(vmult(con.n, con.jnAcc));
	}
	
	return this.swappedColl ? sum : sum.neg();
};

/// Calculate the total impulse including the friction that was applied by this arbiter.
/// This function should only be called from a post-solve, post-step or cpBodyEachArbiter callback.
Arbiter.prototype.totalImpulseWithFriction = function()
{
	var contacts = this.contacts;
	var sum = new Vect(0,0);
	
	for(var i=0, count=contacts.length; i<count; i++){
		var con = contacts[i];
		sum.add(new Vect(con.jnAcc, con.jtAcc).rotate(con.n));
	}

	return this.swappedColl ? sum : sum.neg();
};

/// Calculate the amount of energy lost in a collision including static, but not dynamic friction.
/// This function should only be called from a post-solve, post-step or cpBodyEachArbiter callback.
Arbiter.prototype.totalKE = function()
{
	var eCoef = (1 - this.e)/(1 + this.e);
	var sum = 0;
	
	var contacts = this.contacts;
	for(var i=0, count=contacts.length; i<count; i++){
		var con = contacts[i];
		var jnAcc = con.jnAcc;
		var jtAcc = con.jtAcc;
		
		sum += eCoef*jnAcc*jnAcc/con.nMass + jtAcc*jtAcc/con.tMass;
	}
	
	return sum;
};

/// Causes a collision pair to be ignored as if you returned false from a begin callback.
/// If called from a pre-step callback, you will still need to return false
/// if you want it to be ignored in the current step.
Arbiter.prototype.ignore = function()
{
	this.state = 'ignore';
};

/// Return the colliding shapes involved for this arbiter.
/// The order of their cpSpace.collision_type values will match
/// the order set when the collision handler was registered.
Arbiter.prototype.getA = function()
{
	return this.swappedColl ? this.b : this.a;
};

Arbiter.prototype.getB = function()
{
	return this.swappedColl ? this.a : this.b;
};

/// Returns true if this is the first step a pair of objects started colliding.
Arbiter.prototype.isFirstContact = function()
{
	return this.state === 'first coll';
};

/// A struct that wraps up the important collision data for an arbiter.
var ContactPoint = function(point, normal, dist)
{
	this.point = point;
	this.normal = normal;
	this.dist = dist;
};

/// Return a contact set from an arbiter.
Arbiter.prototype.getContactPointSet = function()
{
	var set = new Array(this.contacts.length);
	
	var i;
	for(i=0; i<set.length; i++){
		set[i] = new ContactPoint(this.contacts[i].p, this.contacts[i].n, this.contacts[i].dist);
	}
	
	return set;
};

/// Get the normal of the @c ith contact point.
Arbiter.prototype.getNormal = function(i)
{
	var n = this.contacts[i].n;
	return this.swappedColl ? vneg(n) : n;
};

/// Get the position of the @c ith contact point.
Arbiter.prototype.getPoint = function(i)
{
	return this.contacts[i].p;
};

/// Get the depth of the @c ith contact point.
Arbiter.prototype.getDepth = function(i)
{
	return this.contacts[i].dist;
};

/*
Arbiter.prototype.threadForBody = function(body)
{
	return (this.body_a === body ? this.thread_a : this.thread_b);
};*/

var unthreadHelper = function(arb, body, prev, next)
{
	// thread_x_y is quite ugly, but it avoids making unnecessary js objects per arbiter.
	if(prev){
		// cpArbiterThreadForBody(prev, body)->next = next;
		if(prev.body_a === body) {
			prev.thread_a_next = next;
		} else {
			prev.thread_b_next = next;
		}
	} else {
		body.arbiterList = next;
	}
	
	if(next){
		// cpArbiterThreadForBody(next, body)->prev = prev;
		if(next.body_a === body){
			next.thread_a_prev = prev;
		} else {
			next.thread_b_prev = prev;
		}
	}
};

Arbiter.prototype.unthread = function()
{
	unthreadHelper(this, this.body_a, this.thread_a_prev, this.thread_a_next);
	unthreadHelper(this, this.body_b, this.thread_b_prev, this.thread_b_next);
	this.thread_a_prev = this.thread_a_next = null;
	this.thread_b_prev = this.thread_b_next = null;
};

//cpFloat
//cpContactsEstimateCrushingImpulse(cpContact *contacts, int numContacts)
//{
//	cpFloat fsum = 0;
//	cpVect vsum = vzero;
//	
//	for(int i=0; i<numContacts; i++){
//		cpContact *con = &contacts[i];
//		cpVect j = vrotate(con.n, v(con.jnAcc, con.jtAcc));
//		
//		fsum += vlength(j);
//		vsum = vadd(vsum, j);
//	}
//	
//	cpFloat vmag = vlength(vsum);
//	return (1 - vmag/fsum);
//}

Arbiter.prototype.update = function(contacts, handler, a, b)
{
	// Arbiters without contact data may exist if a collision function rejected the collision.
	if(this.contacts){
		// Iterate over the possible pairs to look for hash value matches.
		for(var i=0; i<this.contacts.length; i++){
			var old = this.contacts[i];
			
			for(var j=0; j<contacts.length; j++){
				var new_contact = contacts[j];
				
				// This could trigger false positives, but is fairly unlikely nor serious if it does.
				if(new_contact.hash === old.hash){
					// Copy the persistant contact information.
					new_contact.jnAcc = old.jnAcc;
					new_contact.jtAcc = old.jtAcc;
				}
			}
		}
	}
	
	this.contacts = contacts;
	
	this.handler = handler;
	this.swappedColl = (a.collision_type !== handler.a);
	
	this.e = a.e * b.e;
	this.u = a.u * b.u;
	this.surface_vr = vsub(a.surface_v, b.surface_v);
	
	// For collisions between two similar primitive types, the order could have been swapped.
	this.a = a; this.body_a = a.body;
	this.b = b; this.body_b = b.body;
	
	// mark it as new if it's been cached
	if(this.state == 'cached') this.state = 'first coll';
};

Arbiter.prototype.preStep = function(dt, slop, bias)
{
	var a = this.body_a;
	var b = this.body_b;
	
	for(var i=0; i<this.contacts.length; i++){
		var con = this.contacts[i];
		
		// Calculate the offsets.
		con.r1 = vsub(con.p, a.p);
		con.r2 = vsub(con.p, b.p);
		
		// Calculate the mass normal and mass tangent.
		con.nMass = 1/k_scalar(a, b, con.r1, con.r2, con.n);
		con.tMass = 1/k_scalar(a, b, con.r1, con.r2, vperp(con.n));
	
		// Calculate the target bias velocity.
		con.bias = -bias*min(0, con.dist + slop)/dt;
		con.jBias = 0;
		
		// Calculate the target bounce velocity.
		con.bounce = normal_relative_velocity(a, b, con.r1, con.r2, con.n)*this.e;
	}
};

Arbiter.prototype.applyCachedImpulse = function(dt_coef)
{
	if(this.isFirstContact()) return;
	
	var a = this.body_a;
	var b = this.body_b;
	
	for(var i=0; i<this.contacts.length; i++){
		var con = this.contacts[i];
		//var j = vrotate(con.n, new Vect(con.jnAcc, con.jtAcc));
		var nx = con.n.x;
		var ny = con.n.y;
		var jx = nx*con.jnAcc - ny*con.jtAcc;
		var jy = nx*con.jtAcc + ny*con.jnAcc;
		//apply_impulses(a, b, con.r1, con.r2, vmult(j, dt_coef));
		apply_impulses(a, b, con.r1, con.r2, jx * dt_coef, jy * dt_coef);
	}
};

// TODO is it worth splitting velocity/position correction?

var numApplyImpulse = 0;
var numApplyContact = 0;

Arbiter.prototype.applyImpulse = function()
{
	numApplyImpulse++;
	//if (!this.contacts) { throw new Error('contacts is undefined'); }
	var a = this.body_a;
	var b = this.body_b;
	var surface_vr = this.surface_vr;
	var friction = this.u;

	for(var i=0; i<this.contacts.length; i++){
		numApplyContact++;
		var con = this.contacts[i];
		var nMass = con.nMass;
		var n = con.n;
		var r1 = con.r1;
		var r2 = con.r2;
		
		//var vr = relative_velocity(a, b, r1, r2);
		var vrx = b.vx - r2.y * b.w - (a.vx - r1.y * a.w);
		var vry = b.vy + r2.x * b.w - (a.vy + r1.x * a.w);
		
		//var vb1 = vadd(vmult(vperp(r1), a.w_bias), a.v_bias);
		//var vb2 = vadd(vmult(vperp(r2), b.w_bias), b.v_bias);
		//var vbn = vdot(vsub(vb2, vb1), n);

		var vbn = n.x*(b.v_biasx - r2.y * b.w_bias - a.v_biasx + r1.y * a.w_bias) +
				n.y*(r2.x*b.w_bias + b.v_biasy - r1.x * a.w_bias - a.v_biasy);

		var vrn = vdot2(vrx, vry, n.x, n.y);
		//var vrt = vdot(vadd(vr, surface_vr), vperp(n));
		var vrt = vdot2(vrx + surface_vr.x, vry + surface_vr.y, -n.y, n.x);
		
		var jbn = (con.bias - vbn)*nMass;
		var jbnOld = con.jBias;
		con.jBias = max(jbnOld + jbn, 0);
		
		var jn = -(con.bounce + vrn)*nMass;
		var jnOld = con.jnAcc;
		con.jnAcc = max(jnOld + jn, 0);
		
		var jtMax = friction*con.jnAcc;
		var jt = -vrt*con.tMass;
		var jtOld = con.jtAcc;
		con.jtAcc = clamp(jtOld + jt, -jtMax, jtMax);
		
		//apply_bias_impulses(a, b, r1, r2, vmult(n, con.jBias - jbnOld));
		var bias_x = n.x * (con.jBias - jbnOld);
		var bias_y = n.y * (con.jBias - jbnOld);
		apply_bias_impulse(a, -bias_x, -bias_y, r1);
		apply_bias_impulse(b, bias_x, bias_y, r2);

		//apply_impulses(a, b, r1, r2, vrotate(n, new Vect(con.jnAcc - jnOld, con.jtAcc - jtOld)));
		var rot_x = con.jnAcc - jnOld;
		var rot_y = con.jtAcc - jtOld;

		// Inlining apply_impulses decreases speed for some reason :/
		apply_impulses(a, b, r1, r2, n.x*rot_x - n.y*rot_y, n.x*rot_y + n.y*rot_x);
	}
};

Arbiter.prototype.callSeparate = function(space)
{
	// The handler needs to be looked up again as the handler cached on the arbiter may have been deleted since the last step.
	var handler = space.lookupHandler(this.a.collision_type, this.b.collision_type);
	handler.separate(this, space);
};

// From chipmunk_private.h
Arbiter.prototype.next = function(body)
{
	return (this.body_a == body ? this.thread_a_next : this.thread_b_next);
};
/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var numContacts = 0;

var Contact = function(p, n, dist, hash)
{
	this.p = p;
	this.n = n;
	this.dist = dist;
	
	this.r1 = this.r2 = vzero;
	this.nMass = this.tMass = this.bounce = this.bias = 0;

	this.jnAcc = this.jtAcc = this.jBias = 0;
	
	this.hash = hash;
	numContacts++;
};

var NONE = [];

// Add contact points for circle to circle collisions.
// Used by several collision tests.
var circle2circleQuery = function(p1, p2, r1, r2)
{
	var mindist = r1 + r2;
	var delta = vsub(p2, p1);
	var distsq = vlengthsq(delta);
	if(distsq >= mindist*mindist) return;
	
	var dist = Math.sqrt(distsq);

	// Allocate and initialize the contact.
	return new Contact(
		vadd(p1, vmult(delta, 0.5 + (r1 - 0.5*mindist)/(dist ? dist : Infinity))),
		(dist ? vmult(delta, 1/dist) : new Vect(1, 0)),
		dist - mindist,
		0
	);
};

// Collide circle shapes.
var circle2circle = function(circ1, circ2)
{
	var contact = circle2circleQuery(circ1.tc, circ2.tc, circ1.r, circ2.r);
	return contact ? [contact] : NONE;
};

var circle2segment = function(circleShape, segmentShape)
{
	var seg_a = segmentShape.ta;
	var seg_b = segmentShape.tb;
	var center = circleShape.tc;
	
	var seg_delta = vsub(seg_b, seg_a);
	var closest_t = clamp01(vdot(seg_delta, vsub(center, seg_a))/vlengthsq(seg_delta));
	var closest = vadd(seg_a, vmult(seg_delta, closest_t));
	
	var contact = circle2circleQuery(center, closest, circleShape.r, segmentShape.r);
	if(contact){
		var n = contact.n;
		
		// Reject endcap collisions if tangents are provided.
		return (
			(closest_t === 0 && vdot(n, segmentShape.a_tangent) < 0) ||
			(closest_t === 1 && vdot(n, segmentShape.b_tangent) < 0)
		) ? NONE : [contact];
	} else {
		return NONE;
	}
}

// Find the minimum separating axis for the given poly and axis list.
//
// This function needs to return two values - the index of the min. separating axis and
// the value itself. Short of inlining MSA, returning values through a global like this
// is the fastest implementation.
//
// See: http://jsperf.com/return-two-values-from-function/2
var last_MSA_min = 0;
var findMSA = function(poly, planes)
{
	var min_index = 0;
	var min = poly.valueOnAxis(planes[0].n, planes[0].d);
	if(min > 0) return -1;
	
	for(var i=1; i<planes.length; i++){
		var dist = poly.valueOnAxis(planes[i].n, planes[i].d);
		if(dist > 0) {
			return -1;
		} else if(dist > min){
			min = dist;
			min_index = i;
		}
	}
	
	last_MSA_min = min;
	return min_index;
};

// Add contacts for probably penetrating vertexes.
// This handles the degenerate case where an overlap was detected, but no vertexes fall inside
// the opposing polygon. (like a star of david)
var findVertsFallback = function(poly1, poly2, n, dist)
{
	var arr = [];

	var verts1 = poly1.tVerts;
	for(var i=0; i<verts1.length; i+=2){
		var vx = verts1[i];
		var vy = verts1[i+1];
		if(poly2.containsVertPartial(vx, vy, vneg(n))){
			arr.push(new Contact(new Vect(vx, vy), n, dist, hashPair(poly1.hashid, i)));
		}
	}
	
	var verts2 = poly2.tVerts;
	for(var i=0; i<verts2.length; i+=2){
		var vx = verts2[i];
		var vy = verts2[i+1];
		if(poly1.containsVertPartial(vx, vy, n)){
			arr.push(new Contact(new Vect(vx, vy), n, dist, hashPair(poly2.hashid, i)));
		}
	}
	
	return arr;
};

// Add contacts for penetrating vertexes.
var findVerts = function(poly1, poly2, n, dist)
{
	var arr = [];

	var verts1 = poly1.tVerts;
	for(var i=0; i<verts1.length; i+=2){
		var vx = verts1[i];
		var vy = verts1[i+1];
		if(poly2.containsVert(vx, vy)){
			arr.push(new Contact(new Vect(vx, vy), n, dist, hashPair(poly1.hashid, i>>1)));
		}
	}
	
	var verts2 = poly2.tVerts;
	for(var i=0; i<verts2.length; i+=2){
		var vx = verts2[i];
		var vy = verts2[i+1];
		if(poly1.containsVert(vx, vy)){
			arr.push(new Contact(new Vect(vx, vy), n, dist, hashPair(poly2.hashid, i>>1)));
		}
	}
	
	return (arr.length ? arr : findVertsFallback(poly1, poly2, n, dist));
};

// Collide poly shapes together.
var poly2poly = function(poly1, poly2)
{
	var mini1 = findMSA(poly2, poly1.tPlanes);
	if(mini1 == -1) return NONE;
	var min1 = last_MSA_min;
	
	var mini2 = findMSA(poly1, poly2.tPlanes);
	if(mini2 == -1) return NONE;
	var min2 = last_MSA_min;
	
	// There is overlap, find the penetrating verts
	if(min1 > min2)
		return findVerts(poly1, poly2, poly1.tPlanes[mini1].n, min1);
	else
		return findVerts(poly1, poly2, vneg(poly2.tPlanes[mini2].n), min2);
};

// Like cpPolyValueOnAxis(), but for segments.
var segValueOnAxis = function(seg, n, d)
{
	var a = vdot(n, seg.ta) - seg.r;
	var b = vdot(n, seg.tb) - seg.r;
	return min(a, b) - d;
};

// Identify vertexes that have penetrated the segment.
var findPointsBehindSeg = function(arr, seg, poly, pDist, coef) 
{
	var dta = vcross(seg.tn, seg.ta);
	var dtb = vcross(seg.tn, seg.tb);
	var n = vmult(seg.tn, coef);
	
	var verts = poly.tVerts;
	for(var i=0; i<verts.length; i+=2){
		var vx = verts[i];
		var vy = verts[i+1];
		if(vdot2(vx, vy, n.x, n.y) < vdot(seg.tn, seg.ta)*coef + seg.r){
			var dt = vcross2(seg.tn.x, seg.tn.y, vx, vy);
			if(dta >= dt && dt >= dtb){
				arr.push(new Contact(new Vect(vx, vy), n, pDist, hashPair(poly.hashid, i)));
			}
		}
	}
};

// This one is complicated and gross. Just don't go there...
// TODO: Comment me!
var seg2poly = function(seg, poly)
{
	var arr = [];

	var planes = poly.tPlanes;
	var numVerts = planes.length;
	
	var segD = vdot(seg.tn, seg.ta);
	var minNorm = poly.valueOnAxis(seg.tn, segD) - seg.r;
	var minNeg = poly.valueOnAxis(vneg(seg.tn), -segD) - seg.r;
	if(minNeg > 0 || minNorm > 0) return NONE;
	
	var mini = 0;
	var poly_min = segValueOnAxis(seg, planes[0].n, planes[0].d);
	if(poly_min > 0) return NONE;
	for(var i=0; i<numVerts; i++){
		var dist = segValueOnAxis(seg, planes[i].n, planes[i].d);
		if(dist > 0){
			return NONE;
		} else if(dist > poly_min){
			poly_min = dist;
			mini = i;
		}
	}
	
	var poly_n = vneg(planes[mini].n);
	
	var va = vadd(seg.ta, vmult(poly_n, seg.r));
	var vb = vadd(seg.tb, vmult(poly_n, seg.r));
	if(poly.containsVert(va.x, va.y))
		arr.push(new Contact(va, poly_n, poly_min, hashPair(seg.hashid, 0)));
	if(poly.containsVert(vb.x, vb.y))
		arr.push(new Contact(vb, poly_n, poly_min, hashPair(seg.hashid, 1)));
	
	// Floating point precision problems here.
	// This will have to do for now.
//	poly_min -= cp_collision_slop; // TODO is this needed anymore?
	
	if(minNorm >= poly_min || minNeg >= poly_min) {
		if(minNorm > minNeg)
			findPointsBehindSeg(arr, seg, poly, minNorm, 1);
		else
			findPointsBehindSeg(arr, seg, poly, minNeg, -1);
	}
	
	// If no other collision points are found, try colliding endpoints.
	if(arr.length === 0){
		var mini2 = mini * 2;
		var verts = poly.tVerts;

		var poly_a = new Vect(verts[mini2], verts[mini2+1]);
		
		var con;
		if((con = circle2circleQuery(seg.ta, poly_a, seg.r, 0, arr))) return [con];
		if((con = circle2circleQuery(seg.tb, poly_a, seg.r, 0, arr))) return [con];

		var len = numVerts * 2;
		var poly_b = new Vect(verts[(mini2+2)%len], verts[(mini2+3)%len]);
		if((con = circle2circleQuery(seg.ta, poly_b, seg.r, 0, arr))) return [con];
		if((con = circle2circleQuery(seg.tb, poly_b, seg.r, 0, arr))) return [con];
	}

//	console.log(poly.tVerts, poly.tPlanes);
//	console.log('seg2poly', arr);
	return arr;
};

// This one is less gross, but still gross.
// TODO: Comment me!
var circle2poly = function(circ, poly)
{
	var planes = poly.tPlanes;
	
	var mini = 0;
	var min = vdot(planes[0].n, circ.tc) - planes[0].d - circ.r;
	for(var i=0; i<planes.length; i++){
		var dist = vdot(planes[i].n, circ.tc) - planes[i].d - circ.r;
		if(dist > 0){
			return NONE;
		} else if(dist > min) {
			min = dist;
			mini = i;
		}
	}
	
	var n = planes[mini].n;

	var verts = poly.tVerts;
	var len = verts.length;
	var mini2 = mini<<1;

	//var a = poly.tVerts[mini];
	//var b = poly.tVerts[(mini + 1)%poly.tVerts.length];
	var ax = verts[mini2];
	var ay = verts[mini2+1];
	var bx = verts[(mini2+2)%len];
	var by = verts[(mini2+3)%len];

	var dta = vcross2(n.x, n.y, ax, ay);
	var dtb = vcross2(n.x, n.y, bx, by);
	var dt = vcross(n, circ.tc);
		
	if(dt < dtb){
		var con = circle2circleQuery(circ.tc, new Vect(bx, by), circ.r, 0, con);
		return con ? [con] : NONE;
	} else if(dt < dta) {
		return [new Contact(
			vsub(circ.tc, vmult(n, circ.r + min/2)),
			vneg(n),
			min,
			0
		)];
	} else {
		var con = circle2circleQuery(circ.tc, new Vect(ax, ay), circ.r, 0, con);
		return con ? [con] : NONE;
	}
};

// The javascripty way to do this would be either nested object or methods on the prototypes.
// 
// However, the *fastest* way is the method below.
// See: http://jsperf.com/dispatch

// These are copied from the prototypes into the actual objects in the Shape constructor.
CircleShape.prototype.collisionCode = 0;
SegmentShape.prototype.collisionCode = 1;
PolyShape.prototype.collisionCode = 2;

CircleShape.prototype.collisionTable = [
	circle2circle,
	circle2segment,
	circle2poly
];

SegmentShape.prototype.collisionTable = [
	null,
	function(segA, segB) { return NONE; }, // seg2seg
	seg2poly
];

PolyShape.prototype.collisionTable = [
	null,
	null,
	poly2poly
];

var collideShapes = cp.collideShapes = function(a, b)
{
	assert(a.collisionCode <= b.collisionCode, 'Collided shapes must be sorted by type');
	return a.collisionTable[b.collisionCode](a, b);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var defaultCollisionHandler = new CollisionHandler();

/// Basic Unit of Simulation in Chipmunk
var Space = cp.Space = function() {
	this.stamp = 0;
	this.curr_dt = 0;

	this.bodies = [];
	this.rousedBodies = [];
	this.sleepingComponents = [];
	
	this.staticShapes = new BBTree(null);
	this.activeShapes = new BBTree(this.staticShapes);
	
	this.arbiters = [];
	this.contactBuffersHead = null;
	this.cachedArbiters = {};
	//this.pooledArbiters = [];
	
	this.constraints = [];
	
	this.locked = 0;
	
	this.collisionHandlers = {};
	this.defaultHandler = defaultCollisionHandler;

	this.postStepCallbacks = [];
	
	/// Number of iterations to use in the impulse solver to solve contacts.
	this.iterations = 10;
	
	/// Gravity to pass to rigid bodies when integrating velocity.
	this.gravity = vzero;
	
	/// Damping rate expressed as the fraction of velocity bodies retain each second.
	/// A value of 0.9 would mean that each body's velocity will drop 10% per second.
	/// The default value is 1.0, meaning no damping is applied.
	/// @note This damping value is different than those of cpDampedSpring and cpDampedRotarySpring.
	this.damping = 1;
	
	/// Speed threshold for a body to be considered idle.
	/// The default value of 0 means to let the space guess a good threshold based on gravity.
	this.idleSpeedThreshold = 0;
	
	/// Time a group of bodies must remain idle in order to fall asleep.
	/// Enabling sleeping also implicitly enables the the contact graph.
	/// The default value of Infinity disables the sleeping algorithm.
	this.sleepTimeThreshold = Infinity;
	
	/// Amount of encouraged penetration between colliding shapes..
	/// Used to reduce oscillating contacts and keep the collision cache warm.
	/// Defaults to 0.1. If you have poor simulation quality,
	/// increase this number as much as possible without allowing visible amounts of overlap.
	this.collisionSlop = 0.1;
	
	/// Determines how fast overlapping shapes are pushed apart.
	/// Expressed as a fraction of the error remaining after each second.
	/// Defaults to pow(1.0 - 0.1, 60.0) meaning that Chipmunk fixes 10% of overlap each frame at 60Hz.
	this.collisionBias = Math.pow(1 - 0.1, 60);
	
	/// Number of frames that contact information should persist.
	/// Defaults to 3. There is probably never a reason to change this value.
	this.collisionPersistence = 3;
	
	/// Rebuild the contact graph during each step. Must be enabled to use the cpBodyEachArbiter() function.
	/// Disabled by default for a small performance boost. Enabled implicitly when the sleeping feature is enabled.
	this.enableContactGraph = false;
	
	/// The designated static body for this space.
	/// You can modify this body, or replace it with your own static body.
	/// By default it points to a statically allocated cpBody in the cpSpace struct.
	this.staticBody = new Body(Infinity, Infinity);
	this.staticBody.nodeIdleTime = Infinity;

	// Cache the collideShapes callback function for the space.
	this.collideShapes = this.makeCollideShapes();
};

Space.prototype.getCurrentTimeStep = function() { return this.curr_dt; };

Space.prototype.setIterations = function(iter) { this.iterations = iter; };

/// returns true from inside a callback and objects cannot be added/removed.
Space.prototype.isLocked = function()
{
	return this.locked;
};

var assertSpaceUnlocked = function(space)
{
	assert(!space.locked, "This addition/removal cannot be done safely during a call to cpSpaceStep() \
 or during a query. Put these calls into a post-step callback.");
};

// **** Collision handler function management

/// Set a collision handler to be used whenever the two shapes with the given collision types collide.
/// You can pass null for any function you don't want to implement.
Space.prototype.addCollisionHandler = function(a, b, begin, preSolve, postSolve, separate)
{
	assertSpaceUnlocked(this);
		
	// Remove any old function so the new one will get added.
	this.removeCollisionHandler(a, b);
	
	var handler = new CollisionHandler();
	handler.a = a;
	handler.b = b;
	if(begin) handler.begin = begin;
	if(preSolve) handler.preSolve = preSolve;
	if(postSolve) handler.postSolve = postSolve;
	if(separate) handler.separate = separate;

	this.collisionHandlers[hashPair(a, b)] = handler;
};

/// Unset a collision handler.
Space.prototype.removeCollisionHandler = function(a, b)
{
	assertSpaceUnlocked(this);
	
	delete this.collisionHandlers[hashPair(a, b)];
};

/// Set a default collision handler for this space.
/// The default collision handler is invoked for each colliding pair of shapes
/// that isn't explicitly handled by a specific collision handler.
/// You can pass null for any function you don't want to implement.
Space.prototype.setDefaultCollisionHandler = function(begin, preSolve, postSolve, separate)
{
	assertSpaceUnlocked(this);

	var handler = new CollisionHandler();
	if(begin) handler.begin = begin;
	if(preSolve) handler.preSolve = preSolve;
	if(postSolve) handler.postSolve = postSolve;
	if(separate) handler.separate = separate;

	this.defaultHandler = handler;
};

Space.prototype.lookupHandler = function(a, b)
{
	return this.collisionHandlers[hashPair(a, b)] || this.defaultHandler;
};

// **** Body, Shape, and Joint Management

/// Add a collision shape to the simulation.
/// If the shape is attached to a static body, it will be added as a static shape.
Space.prototype.addShape = function(shape)
{
	var body = shape.body;
	if(body.isStatic()) return this.addStaticShape(shape);
	
	assert(!shape.space, "This shape is already added to a space and cannot be added to another.");
	assertSpaceUnlocked(this);
	
	body.activate();
	body.addShape(shape);
	
	shape.update(body.p, body.rot);
	this.activeShapes.insert(shape, shape.hashid);
	shape.space = this;
		
	return shape;
};

/// Explicity add a shape as a static shape to the simulation.
Space.prototype.addStaticShape = function(shape)
{
	assert(!shape.space, "This shape is already added to a space and cannot be added to another.");
	assertSpaceUnlocked(this);
	
	var body = shape.body;
	body.addShape(shape);

	shape.update(body.p, body.rot);
	this.staticShapes.insert(shape, shape.hashid);
	shape.space = this;
	
	return shape;
};

/// Add a rigid body to the simulation.
Space.prototype.addBody = function(body)
{
	assert(!body.isStatic(), "Static bodies cannot be added to a space as they are not meant to be simulated.");
	assert(!body.space, "This body is already added to a space and cannot be added to another.");
	assertSpaceUnlocked(this);
	
	this.bodies.push(body);
	body.space = this;
	
	return body;
};

/// Add a constraint to the simulation.
Space.prototype.addConstraint = function(constraint)
{
	assert(!constraint.space, "This shape is already added to a space and cannot be added to another.");
	assertSpaceUnlocked(this);
	
	var a = constraint.a, b = constraint.b;

	a.activate();
	b.activate();
	this.constraints.push(constraint);
	
	// Push onto the heads of the bodies' constraint lists
	constraint.next_a = a.constraintList; a.constraintList = constraint;
	constraint.next_b = b.constraintList; b.constraintList = constraint;
	constraint.space = this;
	
	return constraint;
};

Space.prototype.filterArbiters = function(body, filter)
{
	for (var hash in this.cachedArbiters)
	{
		var arb = this.cachedArbiters[hash];

		// Match on the filter shape, or if it's null the filter body
		if(
			(body === arb.body_a && (filter === arb.a || filter === null)) ||
			(body === arb.body_b && (filter === arb.b || filter === null))
		){
			// Call separate when removing shapes.
			if(filter && arb.state !== 'cached') arb.callSeparate(this);
			
			arb.unthread();

			deleteObjFromList(this.arbiters, arb);
			//this.pooledArbiters.push(arb);
			
			delete this.cachedArbiters[hash];
		}
	}
};

/// Remove a collision shape from the simulation.
Space.prototype.removeShape = function(shape)
{
	var body = shape.body;
	if(body.isStatic()){
		this.removeStaticShape(shape);
	} else {
		assert(this.containsShape(shape),
			"Cannot remove a shape that was not added to the space. (Removed twice maybe?)");
		assertSpaceUnlocked(this);
		
		body.activate();
		body.removeShape(shape);
		this.filterArbiters(body, shape);
		this.activeShapes.remove(shape, shape.hashid);
		shape.space = null;
	}
};

/// Remove a collision shape added using addStaticShape() from the simulation.
Space.prototype.removeStaticShape = function(shape)
{
	assert(this.containsShape(shape),
		"Cannot remove a static or sleeping shape that was not added to the space. (Removed twice maybe?)");
	assertSpaceUnlocked(this);
	
	var body = shape.body;
	if(body.isStatic()) body.activateStatic(shape);
	body.removeShape(shape);
	this.filterArbiters(body, shape);
	this.staticShapes.remove(shape, shape.hashid);
	shape.space = null;
};

/// Remove a rigid body from the simulation.
Space.prototype.removeBody = function(body)
{
	assert(this.containsBody(body),
		"Cannot remove a body that was not added to the space. (Removed twice maybe?)");
	assertSpaceUnlocked(this);
	
	body.activate();
//	this.filterArbiters(body, null);
	deleteObjFromList(this.bodies, body);
	body.space = null;
};

/// Remove a constraint from the simulation.
Space.prototype.removeConstraint = function(constraint)
{
	assert(this.containsConstraint(constraint),
		"Cannot remove a constraint that was not added to the space. (Removed twice maybe?)");
	assertSpaceUnlocked(this);
	
	constraint.a.activate();
	constraint.b.activate();
	deleteObjFromList(this.constraints, constraint);
	
	constraint.a.removeConstraint(constraint);
	constraint.b.removeConstraint(constraint);
	constraint.space = null;
};

/// Test if a collision shape has been added to the space.
Space.prototype.containsShape = function(shape)
{
	return (shape.space === this);
};

/// Test if a rigid body has been added to the space.
Space.prototype.containsBody = function(body)
{
	return (body.space == this);
};

/// Test if a constraint has been added to the space.
Space.prototype.containsConstraint = function(constraint)
{
	return (constraint.space == this);
};

Space.prototype.uncacheArbiter = function(arb)
{
	delete this.cachedArbiters[hashPair(arb.a.hashid, arb.b.hashid)];
	deleteObjFromList(this.arbiters, arb);
};


// **** Iteration

/// Call @c func for each body in the space.
Space.prototype.eachBody = function(func)
{
	this.lock(); {
		var bodies = this.bodies;
		
		for(var i=0; i<bodies.length; i++){
			func(bodies[i]);
		}
		
		var components = this.sleepingComponents;
		for(var i=0; i<components.length; i++){
			var root = components[i];
			
			var body = root;
			while(body){
				var next = body.nodeNext;
				func(body);
				body = next;
			}
		}
	} this.unlock(true);
};

/// Call @c func for each shape in the space.
Space.prototype.eachShape = function(func)
{
	this.lock(); {
		this.activeShapes.each(func);
		this.staticShapes.each(func);
	} this.unlock(true);
};

/// Call @c func for each shape in the space.
Space.prototype.eachConstraint = function(func)
{
	this.lock(); {
		var constraints = this.constraints;
		
		for(var i=0; i<constraints.length; i++){
			func(constraints[i]);
		}
	} this.unlock(true);
};

// **** Spatial Index Management

/// Update the collision detection info for the static shapes in the space.
Space.prototype.reindexStatic = function()
{
	assert(!this.locked, "You cannot manually reindex objects while the space is locked. Wait until the current query or step is complete.");
	
	this.staticShapes.each(function(shape){
		var body = shape.body;
		shape.update(body.p, body.rot);
	});
	this.staticShapes.reindex();
};

/// Update the collision detection data for a specific shape in the space.
Space.prototype.reindexShape = function(shape)
{
	assert(!this.locked, "You cannot manually reindex objects while the space is locked. Wait until the current query or step is complete.");
	
	var body = shape.body;
	shape.update(body.p, body.rot);
	
	// attempt to rehash the shape in both hashes
	this.activeShapes.reindexObject(shape, shape.hashid);
	this.staticShapes.reindexObject(shape, shape.hashid);
};

/// Update the collision detection data for all shapes attached to a body.
Space.prototype.reindexShapesForBody = function(body)
{
	for(var shape = body.shapeList; shape; shape = shape.next){
		this.reindexShape(shape);
	}
};

/// Switch the space to use a spatial has as it's spatial index.
Space.prototype.useSpatialHash = function(dim, count)
{
	throw new Error('Spatial Hash not yet implemented!');
	
	var staticShapes = new SpaceHash(dim, count, null);
	var activeShapes = new SpaceHash(dim, count, staticShapes);
	
	this.staticShapes.each(function(shape){
		staticShapes.insert(shape, shape.hashid);
	});
	this.activeShapes.each(function(shape){
		activeShapes.insert(shape, shape.hashid);
	});
		
	this.staticShapes = staticShapes;
	this.activeShapes = activeShapes;
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
 
/// **** Sleeping Functions

Space.prototype.activateBody = function(body)
{
	assert(!body.isRogue(), "Internal error: Attempting to activate a rogue body.");
	
	if(this.locked){
		// cpSpaceActivateBody() is called again once the space is unlocked
		if(this.rousedBodies.indexOf(body) === -1) this.rousedBodies.push(body);
	} else {
		this.bodies.push(body);

		for(var i = 0; i < body.shapeList.length; i++){
			var shape = body.shapeList[i];
			this.staticShapes.remove(shape, shape.hashid);
			this.activeShapes.insert(shape, shape.hashid);
		}
		
		for(var arb = body.arbiterList; arb; arb = arb.next(body)){
			var bodyA = arb.body_a;
			if(body === bodyA || bodyA.isStatic()){
				//var contacts = arb.contacts;
				
				// Restore contact values back to the space's contact buffer memory
				//arb.contacts = cpContactBufferGetArray(this);
				//memcpy(arb.contacts, contacts, numContacts*sizeof(cpContact));
				//cpSpacePushContacts(this, numContacts);
				
				// Reinsert the arbiter into the arbiter cache
				var a = arb.a, b = arb.b;
				this.cachedArbiters[hashPair(a.hashid, b.hashid)] = arb;
				
				// Update the arbiter's state
				arb.stamp = this.stamp;
				arb.handler = this.lookupHandler(a.collision_type, b.collision_type);
				this.arbiters.push(arb);
			}
		}
		
		for(var constraint = body.constraintList; constraint; constraint = constraint.nodeNext){
			var bodyA = constraint.a;
			if(body === bodyA || bodyA.isStatic()) this.constraints.push(constraint);
		}
	}
};

Space.prototype.deactivateBody = function(body)
{
	assert(!body.isRogue(), "Internal error: Attempting to deactivate a rogue body.");
	
	deleteObjFromList(this.bodies, body);
	
	for(var i = 0; i < body.shapeList.length; i++){
		var shape = body.shapeList[i];
		this.activeShapes.remove(shape, shape.hashid);
		this.staticShapes.insert(shape, shape.hashid);
	}
	
	for(var arb = body.arbiterList; arb; arb = arb.next(body)){
		var bodyA = arb.body_a;
		if(body === bodyA || bodyA.isStatic()){
			this.uncacheArbiter(arb);
			
			// Save contact values to a new block of memory so they won't time out
			//size_t bytes = arb.numContacts*sizeof(cpContact);
			//cpContact *contacts = (cpContact *)cpcalloc(1, bytes);
			//memcpy(contacts, arb.contacts, bytes);
			//arb.contacts = contacts;
		}
	}
		
	for(var constraint = body.constraintList; constraint; constraint = constraint.nodeNext){
		var bodyA = constraint.a;
		if(body === bodyA || bodyA.isStatic()) deleteObjFromList(this.constraints, constraint);
	}
};

var componentRoot = function(body)
{
	return (body ? body.nodeRoot : null);
};

var componentActivate = function(root)
{
	if(!root || !root.isSleeping(root)) return;
	assert(!root.isRogue(), "Internal Error: componentActivate() called on a rogue body.");
	
	var space = root.space;
	var body = root;
	while(body){
		var next = body.nodeNext;
		
		body.nodeIdleTime = 0;
		body.nodeRoot = null;
		body.nodeNext = null;
		space.activateBody(body);
		
		body = next;
	}
	
	deleteObjFromList(space.sleepingComponents, root);
};

Body.prototype.activate = function()
{
	if(!this.isRogue()){
		this.nodeIdleTime = 0;
		componentActivate(componentRoot(this));
	}
};

Body.prototype.activateStatic = function(filter)
{
	assert(this.isStatic(), "Body.activateStatic() called on a non-static body.");
	
	for(var arb = this.arbiterList; arb; arb = arb.next(this)){
		if(!filter || filter == arb.a || filter == arb.b){
			(arb.body_a == this ? arb.body_b : arb.body_a).activate();
		}
	}
	
	// TODO should also activate joints!
};

Body.prototype.pushArbiter = function(arb)
{
	assertSoft((arb.body_a === this ? arb.thread_a_next : arb.thread_b_next) === null,
		"Internal Error: Dangling contact graph pointers detected. (A)");
	assertSoft((arb.body_a === this ? arb.thread_a_prev : arb.thread_b_prev) === null,
		"Internal Error: Dangling contact graph pointers detected. (B)");
	
	var next = this.arbiterList;
	assertSoft(next === null || (next.body_a === this ? next.thread_a_prev : next.thread_b_prev) === null,
		"Internal Error: Dangling contact graph pointers detected. (C)");

	if(arb.body_a === this){
		arb.thread_a_next = next;
	} else {
		arb.thread_b_next = next;
	}

	if(next){
		if (next.body_a === this){
			next.thread_a_prev = arb;
		} else {
			next.thread_b_prev = arb;
		}
	}
	this.arbiterList = arb;
};

var componentAdd = function(root, body){
	body.nodeRoot = root;

	if(body !== root){
		body.nodeNext = root.nodeNext;
		root.nodeNext = body;
	}
};

var floodFillComponent = function(root, body)
{
	// Rogue bodies cannot be put to sleep and prevent bodies they are touching from sleeping anyway.
	// Static bodies (which are a type of rogue body) are effectively sleeping all the time.
	if(!body.isRogue()){
		var other_root = componentRoot(body);
		if(other_root == null){
			componentAdd(root, body);
			for(var arb = body.arbiterList; arb; arb = arb.next(body)){
				floodFillComponent(root, (body == arb.body_a ? arb.body_b : arb.body_a));
			}
			for(var constraint = body.constraintList; constraint; constraint = constraint.next(body)){
				floodFillComponent(root, (body == constraint.a ? constraint.b : constraint.a));
			}
		} else {
			assertSoft(other_root === root, "Internal Error: Inconsistency detected in the contact graph.");
		}
	}
};

var componentActive = function(root, threshold)
{
	for(var body = root; body; body = body.nodeNext){
		if(body.nodeIdleTime < threshold) return true;
	}
	
	return false;
};

Space.prototype.processComponents = function(dt)
{
	var sleep = (this.sleepTimeThreshold !== Infinity);
	var bodies = this.bodies;

	// These checks can be removed at some stage (if DEBUG == undefined)
	for(var i=0; i<bodies.length; i++){
		var body = bodies[i];
		
		assertSoft(body.nodeNext === null, "Internal Error: Dangling next pointer detected in contact graph.");
		assertSoft(body.nodeRoot === null, "Internal Error: Dangling root pointer detected in contact graph.");
	}

	// Calculate the kinetic energy of all the bodies
	if(sleep){
		var dv = this.idleSpeedThreshold;
		var dvsq = (dv ? dv*dv : vlengthsq(this.gravity)*dt*dt);
	
		for(var i=0; i<bodies.length; i++){
			var body = bodies[i];

			// Need to deal with infinite mass objects
			var keThreshold = (dvsq ? body.m*dvsq : 0);
			body.nodeIdleTime = (body.kineticEnergy() > keThreshold ? 0 : body.nodeIdleTime + dt);
		}
	}

	// Awaken any sleeping bodies found and then push arbiters to the bodies' lists.
	var arbiters = this.arbiters;
	for(var i=0, count=arbiters.length; i<count; i++){
		var arb = arbiters[i];
		var a = arb.body_a, b = arb.body_b;
	
		if(sleep){	
			if((b.isRogue() && !b.isStatic()) || a.isSleeping()) a.activate();
			if((a.isRogue() && !a.isStatic()) || b.isSleeping()) b.activate();
		}
		
		a.pushArbiter(arb);
		b.pushArbiter(arb);
	}
	
	if(sleep){
		// Bodies should be held active if connected by a joint to a non-static rouge body.
		var constraints = this.constraints;
		for(var i=0; i<constraints.length; i++){
			var constraint = constraints[i];
			var a = constraint.a, b = constraint.b;
			
			if(b.isRogue() && !b.isStatic()) a.activate();
			if(a.isRogue() && !a.isStatic()) b.activate();
		}
		
		// Generate components and deactivate sleeping ones
		for(var i=0; i<bodies.length;){
			var body = bodies[i];
			
			if(componentRoot(body) === null){
				// Body not in a component yet. Perform a DFS to flood fill mark 
				// the component in the contact graph using this body as the root.
				floodFillComponent(body, body);
				
				// Check if the component should be put to sleep.
				if(!componentActive(body, this.sleepTimeThreshold)){
					this.sleepingComponents.push(body);
					for(var other = body; other; other = other.nodeNext){
						this.deactivateBody(other);
					}
					
					// deactivateBody() removed the current body from the list.
					// Skip incrementing the index counter.
					continue;
				}
			}
			
			i++;
			
			// Only sleeping bodies retain their component node pointers.
			body.nodeRoot = null;
			body.nodeNext = null;
		}
	}
};

Body.prototype.sleep = function()
{
	this.sleepWithGroup(null);
};

Body.prototype.sleepWithGroup = function(group){
	assert(!this.isStatic() && !this.isRogue(), "Rogue and static bodies cannot be put to sleep.");
	
	var space = this.space;
	assert(space, "Cannot put a rogue body to sleep.");
	assert(!space.locked, "Bodies cannot be put to sleep during a query or a call to cpSpaceStep(). Put these calls into a post-step callback.");
	assert(group === null || group.isSleeping(), "Cannot use a non-sleeping body as a group identifier.");
	
	if(this.isSleeping()){
		assert(componentRoot(this) === componentRoot(group), "The body is already sleeping and it's group cannot be reassigned.");
		return;
	}
	
	for(var i = 0; i < this.shapeList.length; i++){
		this.shapeList[i].update(this.p, this.rot);
	}
	space.deactivateBody(this);
	
	if(group){
		var root = componentRoot(group);
		
		this.nodeRoot = root;
		this.nodeNext = root.nodeNext;
		this.nodeIdleTime = 0;
		
		root.nodeNext = this;
	} else {
		this.nodeRoot = this;
		this.nodeNext = null;
		this.nodeIdleTime = 0;
		
		space.sleepingComponents.push(this);
	}
	
	deleteObjFromList(space.bodies, this);
};

Space.prototype.activateShapesTouchingShape = function(shape){
	if(this.sleepTimeThreshold !== Infinity){
		this.shapeQuery(shape, function(shape, points) {
			shape.body.activate();
		});
	}
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// Point query functions

/// Query the space at a point and call @c func for each shape found.
Space.prototype.pointQuery = function(point, layers, group, func)
{
	var helper = function(shape){
		if(
			!(shape.group && group === shape.group) && (layers & shape.layers) &&
			shape.pointQuery(point)
		){
			func(shape);
		}
	};

	var bb = new BB(point.x, point.y, point.x, point.y);
	this.lock(); {
		this.activeShapes.query(bb, helper);
		this.staticShapes.query(bb, helper);
	} this.unlock(true);
};

/// Query the space at a point and return the first shape found. Returns null if no shapes were found.
Space.prototype.pointQueryFirst = function(point, layers, group)
{
	var outShape = null;
	this.pointQuery(point, layers, group, function(shape) {
		if(!shape.sensor) outShape = shape;
	});
	
	return outShape;
};

// Nearest point query functions

Space.prototype.nearestPointQuery = function(point, maxDistance, layers, group, func)
{
	var helper = function(shape){
		if(!(shape.group && group === shape.group) && (layers & shape.layers)){
			var info = shape.nearestPointQuery(point);

			if(info.d < maxDistance) func(shape, info.d, info.p);
		}
	};

	var bb = bbNewForCircle(point, maxDistance);

	this.lock(); {
		this.activeShapes.query(bb, helper);
		this.staticShapes.query(bb, helper);
	} this.unlock(true);
};

// Unlike the version in chipmunk, this returns a NearestPointQueryInfo object. Use its .shape
// property to get the actual shape.
Space.prototype.nearestPointQueryNearest = function(point, maxDistance, layers, group)
{
	var out;

	var helper = function(shape){
		if(!(shape.group && group === shape.group) && (layers & shape.layers) && !shape.sensor){
			var info = shape.nearestPointQuery(point);

			if(info.d < maxDistance && (!out || info.d < out.d)) out = info;
		}
	};

	var bb = bbNewForCircle(point, maxDistance);
	this.activeShapes.query(bb, helper);
	this.staticShapes.query(bb, helper);

	return out;
};

/// Perform a directed line segment query (like a raycast) against the space calling @c func for each shape intersected.
Space.prototype.segmentQuery = function(start, end, layers, group, func)
{
	var helper = function(shape){
		var info;
		
		if(
			!(shape.group && group === shape.group) && (layers & shape.layers) &&
			(info = shape.segmentQuery(start, end))
		){
			func(shape, info.t, info.n);
		}
		
		return 1;
	};

	this.lock(); {
		this.staticShapes.segmentQuery(start, end, 1, helper);
		this.activeShapes.segmentQuery(start, end, 1, helper);
	} this.unlock(true);
};

/// Perform a directed line segment query (like a raycast) against the space and return the first shape hit.
/// Returns null if no shapes were hit.
Space.prototype.segmentQueryFirst = function(start, end, layers, group)
{
	var out = null;

	var helper = function(shape){
		var info;
		
		if(
			!(shape.group && group === shape.group) && (layers & shape.layers) &&
			!shape.sensor &&
			(info = shape.segmentQuery(start, end)) &&
			(out === null || info.t < out.t)
		){
			out = info;
		}
		
		return out ? out.t : 1;
	};

	this.staticShapes.segmentQuery(start, end, 1, helper);
	this.activeShapes.segmentQuery(start, end, out ? out.t : 1, helper);
	
	return out;
};

/// Perform a fast rectangle query on the space calling @c func for each shape found.
/// Only the shape's bounding boxes are checked for overlap, not their full shape.
Space.prototype.bbQuery = function(bb, layers, group, func)
{
	var helper = function(shape){
		if(
			!(shape.group && group === shape.group) && (layers & shape.layers) &&
			bbIntersects2(bb, shape.bb_l, shape.bb_b, shape.bb_r, shape.bb_t)
		){
			func(shape);
		}
	};
	
	this.lock(); {
		this.activeShapes.query(bb, helper);
		this.staticShapes.query(bb, helper);
	} this.unlock(true);
};

/// Query a space for any shapes overlapping the given shape and call @c func for each shape found.
Space.prototype.shapeQuery = function(shape, func)
{
	var body = shape.body;

	//var bb = (body ? shape.update(body.p, body.rot) : shape.bb);
	if(body){
		shape.update(body.p, body.rot);
	}
	var bb = new BB(shape.bb_l, shape.bb_b, shape.bb_r, shape.bb_t);

	//shapeQueryContext context = {func, data, false};
	var anyCollision = false;
	
	var helper = function(b){
		var a = shape;
		// Reject any of the simple cases
		if(
			(a.group && a.group === b.group) ||
			!(a.layers & b.layers) ||
			a === b
		) return;
		
		var contacts;
		
		// Shape 'a' should have the lower shape type. (required by collideShapes() )
		if(a.collisionCode <= b.collisionCode){
			contacts = collideShapes(a, b);
		} else {
			contacts = collideShapes(b, a);
			for(var i=0; i<contacts.length; i++) contacts[i].n = vneg(contacts[i].n);
		}
		
		if(contacts.length){
			anyCollision = !(a.sensor || b.sensor);
			
			if(func){
				var set = new Array(contacts.length);
				for(var i=0; i<contacts.length; i++){
					set[i] = new ContactPoint(contacts[i].p, contacts[i].n, contacts[i].dist);
				}
				
				func(b, set);
			}
		}
	};

	this.lock(); {
		this.activeShapes.query(bb, helper);
		this.staticShapes.query(bb, helper);
	} this.unlock(true);
	
	return anyCollision;
};

/* Copyright (c) 2007 Scott Lembcke
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// **** Post Step Callback Functions

/// Schedule a post-step callback to be called when cpSpaceStep() finishes.
Space.prototype.addPostStepCallback = function(func)
{
	assertSoft(this.locked,
		"Adding a post-step callback when the space is not locked is unnecessary. " +
		"Post-step callbacks will not called until the end of the next call to cpSpaceStep() or the next query.");

	this.postStepCallbacks.push(func);
};

Space.prototype.runPostStepCallbacks = function()
{
	// Don't cache length because post step callbacks may add more post step callbacks
	// directly or indirectly.
	for(var i = 0; i < this.postStepCallbacks.length; i++){
		this.postStepCallbacks[i]();
	}
	this.postStepCallbacks = [];
};

// **** Locking Functions

Space.prototype.lock = function()
{
	this.locked++;
};

Space.prototype.unlock = function(runPostStep)
{
	this.locked--;
	assert(this.locked >= 0, "Internal Error: Space lock underflow.");

	if(this.locked === 0 && runPostStep){
		var waking = this.rousedBodies;
		for(var i=0; i<waking.length; i++){
			this.activateBody(waking[i]);
		}

		waking.length = 0;

		this.runPostStepCallbacks();
	}
};

// **** Contact Buffer Functions

/* josephg:
 *
 * This code might be faster in JS than just allocating objects each time - I'm
 * really not sure. If the contact buffer solution is used, there will also
 * need to be changes in cpCollision.js to fill a passed array instead of creating
 * new arrays each time.
 *
 * TODO: Benchmark me once chipmunk is working.
 */

/*
var ContactBuffer = function(stamp, splice)
{
	this.stamp = stamp;
	// Contact buffers are a circular linked list.
	this.next = splice ? splice.next : this;
	this.contacts = [];
};

Space.prototype.pushFreshContactBuffer = function()
{
	var stamp = this.stamp;

	var head = this.contactBuffersHead;

	if(!head){
		// No buffers have been allocated, make one
		this.contactBuffersHead = new ContactBuffer(stamp, null);
	} else if(stamp - head.next.stamp > this.collisionPersistence){
		// The tail buffer is available, rotate the ring
		var tail = head.next;
		tail.stamp = stamp;
		tail.contacts.length = 0;
		this.contactBuffersHead = tail;
	} else {
		// Allocate a new buffer and push it into the ring
		var buffer = new ContactBuffer(stamp, head);
		this.contactBuffersHead = head.next = buffer;
	}
};

cpContact *
cpContactBufferGetArray(cpSpace *space)
{
	if(space.contactBuffersHead.numContacts + CP_MAX_CONTACTS_PER_ARBITER > CP_CONTACTS_BUFFER_SIZE){
		// contact buffer could overflow on the next collision, push a fresh one.
		space.pushFreshContactBuffer();
	}

	cpContactBufferHeader *head = space.contactBuffersHead;
	return ((cpContactBuffer *)head)->contacts + head.numContacts;
}

void
cpSpacePushContacts(cpSpace *space, int count)
{
	cpAssertHard(count <= CP_MAX_CONTACTS_PER_ARBITER, "Internal Error: Contact buffer overflow!");
	space.contactBuffersHead.numContacts += count;
}

static void
cpSpacePopContacts(cpSpace *space, int count){
	space.contactBuffersHead.numContacts -= count;
}
*/

// **** Collision Detection Functions

/* Use this to re-enable object pools.
static void *
cpSpaceArbiterSetTrans(cpShape **shapes, cpSpace *space)
{
	if(space.pooledArbiters.num == 0){
		// arbiter pool is exhausted, make more
		int count = CP_BUFFER_BYTES/sizeof(cpArbiter);
		cpAssertHard(count, "Internal Error: Buffer size too small.");

		cpArbiter *buffer = (cpArbiter *)cpcalloc(1, CP_BUFFER_BYTES);
		cpArrayPush(space.allocatedBuffers, buffer);

		for(int i=0; i<count; i++) cpArrayPush(space.pooledArbiters, buffer + i);
	}

	return cpArbiterInit((cpArbiter *)cpArrayPop(space.pooledArbiters), shapes[0], shapes[1]);
}*/

// Callback from the spatial hash.
Space.prototype.makeCollideShapes = function()
{
	// It would be nicer to use .bind() or something, but this is faster.
	var space_ = this;
	return function(a, b){
		var space = space_;

		// Reject any of the simple cases
		if(
			// BBoxes must overlap
			//!bbIntersects(a.bb, b.bb)
			!(a.bb_l <= b.bb_r && b.bb_l <= a.bb_r && a.bb_b <= b.bb_t && b.bb_b <= a.bb_t)
			// Don't collide shapes attached to the same body.
			|| a.body === b.body
			// Don't collide objects in the same non-zero group
			|| (a.group && a.group === b.group)
			// Don't collide objects that don't share at least on layer.
			|| !(a.layers & b.layers)
		) return;

		var handler = space.lookupHandler(a.collision_type, b.collision_type);

		var sensor = a.sensor || b.sensor;
		if(sensor && handler === defaultCollisionHandler) return;

		// Shape 'a' should have the lower shape type. (required by cpCollideShapes() )
		if(a.collisionCode > b.collisionCode){
			var temp = a;
			a = b;
			b = temp;
		}

		// Narrow-phase collision detection.
		//cpContact *contacts = cpContactBufferGetArray(space);
		//int numContacts = cpCollideShapes(a, b, contacts);
		var contacts = collideShapes(a, b);
		if(contacts.length === 0) return; // Shapes are not colliding.
		//cpSpacePushContacts(space, numContacts);

		// Get an arbiter from space.arbiterSet for the two shapes.
		// This is where the persistant contact magic comes from.
		var arbHash = hashPair(a.hashid, b.hashid);
		var arb = space.cachedArbiters[arbHash];
		if (!arb){
			arb = space.cachedArbiters[arbHash] = new Arbiter(a, b);
		}

		arb.update(contacts, handler, a, b);

		// Call the begin function first if it's the first step
		if(arb.state == 'first coll' && !handler.begin(arb, space)){
			arb.ignore(); // permanently ignore the collision until separation
		}

		if(
			// Ignore the arbiter if it has been flagged
			(arb.state !== 'ignore') &&
			// Call preSolve
			handler.preSolve(arb, space) &&
			// Process, but don't add collisions for sensors.
			!sensor
		){
			space.arbiters.push(arb);
		} else {
			//cpSpacePopContacts(space, numContacts);

			arb.contacts = null;

			// Normally arbiters are set as used after calling the post-solve callback.
			// However, post-solve callbacks are not called for sensors or arbiters rejected from pre-solve.
			if(arb.state !== 'ignore') arb.state = 'normal';
		}

		// Time stamp the arbiter so we know it was used recently.
		arb.stamp = space.stamp;
	};
};

// Hashset filter func to throw away old arbiters.
Space.prototype.arbiterSetFilter = function(arb)
{
	var ticks = this.stamp - arb.stamp;

	var a = arb.body_a, b = arb.body_b;

	// TODO should make an arbiter state for this so it doesn't require filtering arbiters for
	// dangling body pointers on body removal.
	// Preserve arbiters on sensors and rejected arbiters for sleeping objects.
	// This prevents errant separate callbacks from happenening.
	if(
		(a.isStatic() || a.isSleeping()) &&
		(b.isStatic() || b.isSleeping())
	){
		return true;
	}

	// Arbiter was used last frame, but not this one
	if(ticks >= 1 && arb.state != 'cached'){
		arb.callSeparate(this);
		arb.state = 'cached';
	}

	if(ticks >= this.collisionPersistence){
		arb.contacts = null;

		//cpArrayPush(this.pooledArbiters, arb);
		return false;
	}

	return true;
};

// **** All Important cpSpaceStep() Function

var updateFunc = function(shape)
{
	var body = shape.body;
	shape.update(body.p, body.rot);
};

/// Step the space forward in time by @c dt.
Space.prototype.step = function(dt)
{
	// don't step if the timestep is 0!
	if(dt === 0) return;

	assert(vzero.x === 0 && vzero.y === 0, "vzero is invalid");

	this.stamp++;

	var prev_dt = this.curr_dt;
	this.curr_dt = dt;

    var i;
    var j;
    var hash;
	var bodies = this.bodies;
	var constraints = this.constraints;
	var arbiters = this.arbiters;

	// Reset and empty the arbiter lists.
	for(i=0; i<arbiters.length; i++){
		var arb = arbiters[i];
		arb.state = 'normal';

		// If both bodies are awake, unthread the arbiter from the contact graph.
		if(!arb.body_a.isSleeping() && !arb.body_b.isSleeping()){
			arb.unthread();
		}
	}
	arbiters.length = 0;

	this.lock(); {
		// Integrate positions
		for(i=0; i<bodies.length; i++){
			bodies[i].position_func(dt);
		}

		// Find colliding pairs.
		//this.pushFreshContactBuffer();
		this.activeShapes.each(updateFunc);
		this.activeShapes.reindexQuery(this.collideShapes);
	} this.unlock(false);

	// Rebuild the contact graph (and detect sleeping components if sleeping is enabled)
	this.processComponents(dt);

	this.lock(); {
		// Clear out old cached arbiters and call separate callbacks
		for(hash in this.cachedArbiters) {
			if(!this.arbiterSetFilter(this.cachedArbiters[hash])) {
				delete this.cachedArbiters[hash];
			}
		}

		// Prestep the arbiters and constraints.
		var slop = this.collisionSlop;
		var biasCoef = 1 - Math.pow(this.collisionBias, dt);
		for(i=0; i<arbiters.length; i++){
			arbiters[i].preStep(dt, slop, biasCoef);
		}

		for(i=0; i<constraints.length; i++){
			var constraint = constraints[i];

			constraint.preSolve(this);
			constraint.preStep(dt);
		}

		// Integrate velocities.
		var damping = Math.pow(this.damping, dt);
		var gravity = this.gravity;
		for(i=0; i<bodies.length; i++){
			bodies[i].velocity_func(gravity, damping, dt);
		}

		// Apply cached impulses
		var dt_coef = (prev_dt === 0 ? 0 : dt/prev_dt);
		for(i=0; i<arbiters.length; i++){
			arbiters[i].applyCachedImpulse(dt_coef);
		}

		for(i=0; i<constraints.length; i++){
			constraints[i].applyCachedImpulse(dt_coef);
		}

		// Run the impulse solver.
		for(i=0; i<this.iterations; i++){
			for(j=0; j<arbiters.length; j++){
				arbiters[j].applyImpulse();
			}

			for(j=0; j<constraints.length; j++){
				constraints[j].applyImpulse();
			}
		}

		// Run the constraint post-solve callbacks
		for(i=0; i<constraints.length; i++){
			constraints[i].postSolve(this);
		}

		// run the post-solve callbacks
		for(i=0; i<arbiters.length; i++){
			arbiters[i].handler.postSolve(arbiters[i], this);
		}
	} this.unlock(true);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// These are utility routines to use when creating custom constraints.
// I'm not sure if this should be part of the private API or not.
// I should probably clean up the naming conventions if it is...

//#define J_MAX(constraint, dt) (((cpConstraint *)constraint)->maxForce*(dt))

// a and b are bodies.
var relative_velocity = function(a, b, r1, r2){
	//var v1_sum = vadd(a.v, vmult(vperp(r1), a.w));
	var v1_sumx = a.vx + (-r1.y) * a.w;
	var v1_sumy = a.vy + ( r1.x) * a.w;

	//var v2_sum = vadd(b.v, vmult(vperp(r2), b.w));
	var v2_sumx = b.vx + (-r2.y) * b.w;
	var v2_sumy = b.vy + ( r2.x) * b.w;
	
//	return vsub(v2_sum, v1_sum);
	return new Vect(v2_sumx - v1_sumx, v2_sumy - v1_sumy);
};

var normal_relative_velocity = function(a, b, r1, r2, n){
	//return vdot(relative_velocity(a, b, r1, r2), n);
	var v1_sumx = a.vx + (-r1.y) * a.w;
	var v1_sumy = a.vy + ( r1.x) * a.w;
	var v2_sumx = b.vx + (-r2.y) * b.w;
	var v2_sumy = b.vy + ( r2.x) * b.w;

	return vdot2(v2_sumx - v1_sumx, v2_sumy - v1_sumy, n.x, n.y);
};

/*
var apply_impulse = function(body, j, r){
	body.v = vadd(body.v, vmult(j, body.m_inv));
	body.w += body.i_inv*vcross(r, j);
};

var apply_impulses = function(a, b, r1, r2, j)
{
	apply_impulse(a, vneg(j), r1);
	apply_impulse(b, j, r2);
};
*/

var apply_impulse = function(body, jx, jy, r){
//	body.v = body.v.add(vmult(j, body.m_inv));
	body.vx += jx * body.m_inv;
	body.vy += jy * body.m_inv;
//	body.w += body.i_inv*vcross(r, j);
	body.w += body.i_inv*(r.x*jy - r.y*jx);
};

var apply_impulses = function(a, b, r1, r2, jx, jy)
{
	apply_impulse(a, -jx, -jy, r1);
	apply_impulse(b, jx, jy, r2);
};

var apply_bias_impulse = function(body, jx, jy, r)
{
	//body.v_bias = vadd(body.v_bias, vmult(j, body.m_inv));
	body.v_biasx += jx * body.m_inv;
	body.v_biasy += jy * body.m_inv;
	body.w_bias += body.i_inv*vcross2(r.x, r.y, jx, jy);
};

/*
var apply_bias_impulses = function(a, b, r1, r2, j)
{
	apply_bias_impulse(a, vneg(j), r1);
	apply_bias_impulse(b, j, r2);
};*/

var k_scalar_body = function(body, r, n)
{
	var rcn = vcross(r, n);
	return body.m_inv + body.i_inv*rcn*rcn;
};

var k_scalar = function(a, b, r1, r2, n)
{
	var value = k_scalar_body(a, r1, n) + k_scalar_body(b, r2, n);
	assertSoft(value !== 0, "Unsolvable collision or constraint.");
	
	return value;
};

// k1 and k2 are modified by the function to contain the outputs.
var k_tensor = function(a, b, r1, r2, k1, k2)
{
	// calculate mass matrix
	// If I wasn't lazy and wrote a proper matrix class, this wouldn't be so gross...
	var k11, k12, k21, k22;
	var m_sum = a.m_inv + b.m_inv;
	
	// start with I*m_sum
	k11 = m_sum; k12 = 0;
	k21 = 0;     k22 = m_sum;
	
	// add the influence from r1
	var a_i_inv = a.i_inv;
	var r1xsq =  r1.x * r1.x * a_i_inv;
	var r1ysq =  r1.y * r1.y * a_i_inv;
	var r1nxy = -r1.x * r1.y * a_i_inv;
	k11 += r1ysq; k12 += r1nxy;
	k21 += r1nxy; k22 += r1xsq;
	
	// add the influnce from r2
	var b_i_inv = b.i_inv;
	var r2xsq =  r2.x * r2.x * b_i_inv;
	var r2ysq =  r2.y * r2.y * b_i_inv;
	var r2nxy = -r2.x * r2.y * b_i_inv;
	k11 += r2ysq; k12 += r2nxy;
	k21 += r2nxy; k22 += r2xsq;
	
	// invert
	var determinant = k11*k22 - k12*k21;
	assertSoft(determinant !== 0, "Unsolvable constraint.");
	
	var det_inv = 1/determinant;

	k1.x =  k22*det_inv; k1.y = -k12*det_inv;
	k2.x = -k21*det_inv; k2.y =  k11*det_inv;
};

var mult_k = function(vr, k1, k2)
{
	return new Vect(vdot(vr, k1), vdot(vr, k2));
};

var bias_coef = function(errorBias, dt)
{
	return 1 - Math.pow(errorBias, dt);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// TODO: Comment me!

// a and b are bodies that the constraint applies to.
var Constraint = cp.Constraint = function(a, b)
{
	/// The first body connected to this constraint.
	this.a = a;
	/// The second body connected to this constraint.
	this.b = b;
	
	this.space = null;

	this.next_a = null;
	this.next_b = null;
	
	/// The maximum force that this constraint is allowed to use.
	this.maxForce = Infinity;
	/// The rate at which joint error is corrected.
	/// Defaults to pow(1 - 0.1, 60) meaning that it will
	/// correct 10% of the error every 1/60th of a second.
	this.errorBias = Math.pow(1 - 0.1, 60);
	/// The maximum rate at which joint error is corrected.
	this.maxBias = Infinity;
};

Constraint.prototype.activateBodies = function()
{
	if(this.a) this.a.activate();
	if(this.b) this.b.activate();
};

/// These methods are overridden by the constraint itself.
Constraint.prototype.preStep = function(dt) {};
Constraint.prototype.applyCachedImpulse = function(dt_coef) {};
Constraint.prototype.applyImpulse = function() {};
Constraint.prototype.getImpulse = function() { return 0; };

/// Function called before the solver runs. This can be overridden by the user
/// to customize the constraint.
/// Animate your joint anchors, update your motor torque, etc.
Constraint.prototype.preSolve = function(space) {};

/// Function called after the solver runs. This can be overridden by the user
/// to customize the constraint.
/// Use the applied impulse to perform effects like breakable joints.
Constraint.prototype.postSolve = function(space) {};

Constraint.prototype.next = function(body)
{
	return (this.a === body ? this.next_a : this.next_b);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var PinJoint = cp.PinJoint = function(a, b, anchr1, anchr2)
{
	Constraint.call(this, a, b);
	
	this.anchr1 = anchr1;
	this.anchr2 = anchr2;
	
	// STATIC_BODY_CHECK
	var p1 = (a ? vadd(a.p, vrotate(anchr1, a.rot)) : anchr1);
	var p2 = (b ? vadd(b.p, vrotate(anchr2, b.rot)) : anchr2);
	this.dist = vlength(vsub(p2, p1));
	
	assertSoft(this.dist > 0, "You created a 0 length pin joint. A pivot joint will be much more stable.");

	this.r1 = this.r2 = null;
	this.n = null;
	this.nMass = 0;

	this.jnAcc = this.jnMax = 0;
	this.bias = 0;
};

PinJoint.prototype = Object.create(Constraint.prototype);

PinJoint.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	this.r1 = vrotate(this.anchr1, a.rot);
	this.r2 = vrotate(this.anchr2, b.rot);
	
	var delta = vsub(vadd(b.p, this.r2), vadd(a.p, this.r1));
	var dist = vlength(delta);
	this.n = vmult(delta, 1/(dist ? dist : Infinity));
	
	// calculate mass normal
	this.nMass = 1/k_scalar(a, b, this.r1, this.r2, this.n);
	
	// calculate bias velocity
	var maxBias = this.maxBias;
	this.bias = clamp(-bias_coef(this.errorBias, dt)*(dist - this.dist)/dt, -maxBias, maxBias);
	
	// compute max impulse
	this.jnMax = this.maxForce * dt;
};

PinJoint.prototype.applyCachedImpulse = function(dt_coef)
{
	var j = vmult(this.n, this.jnAcc*dt_coef);
	apply_impulses(this.a, this.b, this.r1, this.r2, j.x, j.y);
};

PinJoint.prototype.applyImpulse = function()
{
	var a = this.a;
	var b = this.b;
	var n = this.n;

	// compute relative velocity
	var vrn = normal_relative_velocity(a, b, this.r1, this.r2, n);
	
	// compute normal impulse
	var jn = (this.bias - vrn)*this.nMass;
	var jnOld = this.jnAcc;
	this.jnAcc = clamp(jnOld + jn, -this.jnMax, this.jnMax);
	jn = this.jnAcc - jnOld;
	
	// apply impulse
	apply_impulses(a, b, this.r1, this.r2, n.x*jn, n.y*jn);
};

PinJoint.prototype.getImpulse = function()
{
	return Math.abs(this.jnAcc);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var SlideJoint = cp.SlideJoint = function(a, b, anchr1, anchr2, min, max)
{
	Constraint.call(this, a, b);
	
	this.anchr1 = anchr1;
	this.anchr2 = anchr2;
	this.min = min;
	this.max = max;

	this.r1 = this.r2 = this.n = null;
	this.nMass = 0;
	
	this.jnAcc = this.jnMax = 0;
	this.bias = 0;
};

SlideJoint.prototype = Object.create(Constraint.prototype);

SlideJoint.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	this.r1 = vrotate(this.anchr1, a.rot);
	this.r2 = vrotate(this.anchr2, b.rot);
	
	var delta = vsub(vadd(b.p, this.r2), vadd(a.p, this.r1));
	var dist = vlength(delta);
	var pdist = 0;
	if(dist > this.max) {
		pdist = dist - this.max;
		this.n = vnormalize_safe(delta);
	} else if(dist < this.min) {
		pdist = this.min - dist;
		this.n = vneg(vnormalize_safe(delta));
	} else {
		this.n = vzero;
		this.jnAcc = 0;
	}
	
	// calculate mass normal
	this.nMass = 1/k_scalar(a, b, this.r1, this.r2, this.n);
	
	// calculate bias velocity
	var maxBias = this.maxBias;
	this.bias = clamp(-bias_coef(this.errorBias, dt)*pdist/dt, -maxBias, maxBias);
	
	// compute max impulse
	this.jnMax = this.maxForce * dt;
};

SlideJoint.prototype.applyCachedImpulse = function(dt_coef)
{
	var jn = this.jnAcc * dt_coef;
	apply_impulses(this.a, this.b, this.r1, this.r2, this.n.x * jn, this.n.y * jn);
};

SlideJoint.prototype.applyImpulse = function()
{
	if(this.n.x === 0 && this.n.y === 0) return;  // early exit

	var a = this.a;
	var b = this.b;
	
	var n = this.n;
	var r1 = this.r1;
	var r2 = this.r2;
		
	// compute relative velocity
	var vr = relative_velocity(a, b, r1, r2);
	var vrn = vdot(vr, n);
	
	// compute normal impulse
	var jn = (this.bias - vrn)*this.nMass;
	var jnOld = this.jnAcc;
	this.jnAcc = clamp(jnOld + jn, -this.jnMax, 0);
	jn = this.jnAcc - jnOld;
	
	// apply impulse
	apply_impulses(a, b, this.r1, this.r2, n.x * jn, n.y * jn);
};

SlideJoint.prototype.getImpulse = function()
{
	return Math.abs(this.jnAcc);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// Pivot joints can also be created with (a, b, pivot);
var PivotJoint = cp.PivotJoint = function(a, b, anchr1, anchr2)
{
	Constraint.call(this, a, b);
	
	if(typeof anchr2 === 'undefined') {
		var pivot = anchr1;

		anchr1 = (a ? a.world2Local(pivot) : pivot);
		anchr2 = (b ? b.world2Local(pivot) : pivot);
	}

	this.anchr1 = anchr1;
	this.anchr2 = anchr2;

	this.r1 = this.r2 = vzero;
	
	this.k1 = new Vect(0,0); this.k2 = new Vect(0,0);

	this.jAcc = vzero;

	this.jMaxLen = 0;
	this.bias = vzero;
};

PivotJoint.prototype = Object.create(Constraint.prototype);

PivotJoint.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	this.r1 = vrotate(this.anchr1, a.rot);
	this.r2 = vrotate(this.anchr2, b.rot);
	
	// Calculate mass tensor. Result is stored into this.k1 & this.k2.
	k_tensor(a, b, this.r1, this.r2, this.k1, this.k2);
	
	// compute max impulse
	this.jMaxLen = this.maxForce * dt;
	
	// calculate bias velocity
	var delta = vsub(vadd(b.p, this.r2), vadd(a.p, this.r1));
	this.bias = vclamp(vmult(delta, -bias_coef(this.errorBias, dt)/dt), this.maxBias);
};

PivotJoint.prototype.applyCachedImpulse = function(dt_coef)
{
	apply_impulses(this.a, this.b, this.r1, this.r2, this.jAcc.x * dt_coef, this.jAcc.y * dt_coef);
};

PivotJoint.prototype.applyImpulse = function()
{
	var a = this.a;
	var b = this.b;
	
	var r1 = this.r1;
	var r2 = this.r2;
		
	// compute relative velocity
	var vr = relative_velocity(a, b, r1, r2);
	
	// compute normal impulse
	var j = mult_k(vsub(this.bias, vr), this.k1, this.k2);
	var jOld = this.jAcc;
	this.jAcc = vclamp(vadd(this.jAcc, j), this.jMaxLen);
	
	// apply impulse
	apply_impulses(a, b, this.r1, this.r2, this.jAcc.x - jOld.x, this.jAcc.y - jOld.y);
};

PivotJoint.prototype.getImpulse = function()
{
	return vlength(this.jAcc);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var GrooveJoint = cp.GrooveJoint = function(a, b, groove_a, groove_b, anchr2)
{
	Constraint.call(this, a, b);
	
	this.grv_a = groove_a;
	this.grv_b = groove_b;
	this.grv_n = vperp(vnormalize(vsub(groove_b, groove_a)));
	this.anchr2 = anchr2;
	
	this.grv_tn = null;
	this.clamp = 0;
	this.r1 = this.r2 = null;

	this.k1 = new Vect(0,0);
	this.k2 = new Vect(0,0);

	this.jAcc = vzero;
	this.jMaxLen = 0;
	this.bias = null;
};

GrooveJoint.prototype = Object.create(Constraint.prototype);

GrooveJoint.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	// calculate endpoints in worldspace
	var ta = a.local2World(this.grv_a);
	var tb = a.local2World(this.grv_b);

	// calculate axis
	var n = vrotate(this.grv_n, a.rot);
	var d = vdot(ta, n);
	
	this.grv_tn = n;
	this.r2 = vrotate(this.anchr2, b.rot);
	
	// calculate tangential distance along the axis of r2
	var td = vcross(vadd(b.p, this.r2), n);
	// calculate clamping factor and r2
	if(td <= vcross(ta, n)){
		this.clamp = 1;
		this.r1 = vsub(ta, a.p);
	} else if(td >= vcross(tb, n)){
		this.clamp = -1;
		this.r1 = vsub(tb, a.p);
	} else {
		this.clamp = 0;
		this.r1 = vsub(vadd(vmult(vperp(n), -td), vmult(n, d)), a.p);
	}
	
	// Calculate mass tensor
	k_tensor(a, b, this.r1, this.r2, this.k1, this.k2);	
	
	// compute max impulse
	this.jMaxLen = this.maxForce * dt;
	
	// calculate bias velocity
	var delta = vsub(vadd(b.p, this.r2), vadd(a.p, this.r1));
	this.bias = vclamp(vmult(delta, -bias_coef(this.errorBias, dt)/dt), this.maxBias);
};

GrooveJoint.prototype.applyCachedImpulse = function(dt_coef)
{
	apply_impulses(this.a, this.b, this.r1, this.r2, this.jAcc.x * dt_coef, this.jAcc.y * dt_coef);
};

GrooveJoint.prototype.grooveConstrain = function(j){
	var n = this.grv_tn;
	var jClamp = (this.clamp*vcross(j, n) > 0) ? j : vproject(j, n);
	return vclamp(jClamp, this.jMaxLen);
};

GrooveJoint.prototype.applyImpulse = function()
{
	var a = this.a;
	var b = this.b;
	
	var r1 = this.r1;
	var r2 = this.r2;
	
	// compute impulse
	var vr = relative_velocity(a, b, r1, r2);

	var j = mult_k(vsub(this.bias, vr), this.k1, this.k2);
	var jOld = this.jAcc;
	this.jAcc = this.grooveConstrain(vadd(jOld, j));
	
	// apply impulse
	apply_impulses(a, b, this.r1, this.r2, this.jAcc.x - jOld.x, this.jAcc.y - jOld.y);
};

GrooveJoint.prototype.getImpulse = function()
{
	return vlength(this.jAcc);
};

GrooveJoint.prototype.setGrooveA = function(value)
{
	this.grv_a = value;
	this.grv_n = vperp(vnormalize(vsub(this.grv_b, value)));
	
	this.activateBodies();
};

GrooveJoint.prototype.setGrooveB = function(value)
{
	this.grv_b = value;
	this.grv_n = vperp(vnormalize(vsub(value, this.grv_a)));
	
	this.activateBodies();
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var defaultSpringForce = function(spring, dist){
	return (spring.restLength - dist)*spring.stiffness;
};

var DampedSpring = cp.DampedSpring = function(a, b, anchr1, anchr2, restLength, stiffness, damping)
{
	Constraint.call(this, a, b);
	
	this.anchr1 = anchr1;
	this.anchr2 = anchr2;
	
	this.restLength = restLength;
	this.stiffness = stiffness;
	this.damping = damping;
	this.springForceFunc = defaultSpringForce;

	this.target_vrn = this.v_coef = 0;

	this.r1 = this.r2 = null;
	this.nMass = 0;
	this.n = null;
};

DampedSpring.prototype = Object.create(Constraint.prototype);

DampedSpring.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	this.r1 = vrotate(this.anchr1, a.rot);
	this.r2 = vrotate(this.anchr2, b.rot);
	
	var delta = vsub(vadd(b.p, this.r2), vadd(a.p, this.r1));
	var dist = vlength(delta);
	this.n = vmult(delta, 1/(dist ? dist : Infinity));
	
	var k = k_scalar(a, b, this.r1, this.r2, this.n);
	assertSoft(k !== 0, "Unsolvable this.");
	this.nMass = 1/k;
	
	this.target_vrn = 0;
	this.v_coef = 1 - Math.exp(-this.damping*dt*k);

	// apply this force
	var f_spring = this.springForceFunc(this, dist);
	apply_impulses(a, b, this.r1, this.r2, this.n.x * f_spring * dt, this.n.y * f_spring * dt);
};

DampedSpring.prototype.applyCachedImpulse = function(dt_coef){};

DampedSpring.prototype.applyImpulse = function()
{
	var a = this.a;
	var b = this.b;
	
	var n = this.n;
	var r1 = this.r1;
	var r2 = this.r2;

	// compute relative velocity
	var vrn = normal_relative_velocity(a, b, r1, r2, n);
	
	// compute velocity loss from drag
	var v_damp = (this.target_vrn - vrn)*this.v_coef;
	this.target_vrn = vrn + v_damp;
	
	v_damp *= this.nMass;
	apply_impulses(a, b, this.r1, this.r2, this.n.x * v_damp, this.n.y * v_damp);
};

DampedSpring.prototype.getImpulse = function()
{
	return 0;
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var defaultSpringTorque = function(spring, relativeAngle){
	return (relativeAngle - spring.restAngle)*spring.stiffness;
}

var DampedRotarySpring = cp.DampedRotarySpring = function(a, b, restAngle, stiffness, damping)
{
	Constraint.call(this, a, b);
	
	this.restAngle = restAngle;
	this.stiffness = stiffness;
	this.damping = damping;
	this.springTorqueFunc = defaultSpringTorque;

	this.target_wrn = 0;
	this.w_coef = 0;
	this.iSum = 0;
};

DampedRotarySpring.prototype = Object.create(Constraint.prototype);

DampedRotarySpring.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	var moment = a.i_inv + b.i_inv;
	assertSoft(moment !== 0, "Unsolvable spring.");
	this.iSum = 1/moment;

	this.w_coef = 1 - Math.exp(-this.damping*dt*moment);
	this.target_wrn = 0;

	// apply this torque
	var j_spring = this.springTorqueFunc(this, a.a - b.a)*dt;
	a.w -= j_spring*a.i_inv;
	b.w += j_spring*b.i_inv;
};

// DampedRotarySpring.prototype.applyCachedImpulse = function(dt_coef){};

DampedRotarySpring.prototype.applyImpulse = function()
{
	var a = this.a;
	var b = this.b;
	
	// compute relative velocity
	var wrn = a.w - b.w;//normal_relative_velocity(a, b, r1, r2, n) - this.target_vrn;
	
	// compute velocity loss from drag
	// not 100% certain spring is derived correctly, though it makes sense
	var w_damp = (this.target_wrn - wrn)*this.w_coef;
	this.target_wrn = wrn + w_damp;
	
	//apply_impulses(a, b, this.r1, this.r2, vmult(this.n, v_damp*this.nMass));
	var j_damp = w_damp*this.iSum;
	a.w += j_damp*a.i_inv;
	b.w -= j_damp*b.i_inv;
};

// DampedRotarySpring.prototype.getImpulse = function(){ return 0; };

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var RotaryLimitJoint = cp.RotaryLimitJoint = function(a, b, min, max)
{
	Constraint.call(this, a, b);
	
	this.min = min;
	this.max = max;

	this.jAcc = 0;

	this.iSum = this.bias = this.jMax = 0;
};

RotaryLimitJoint.prototype = Object.create(Constraint.prototype);

RotaryLimitJoint.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	var dist = b.a - a.a;
	var pdist = 0;
	if(dist > this.max) {
		pdist = this.max - dist;
	} else if(dist < this.min) {
		pdist = this.min - dist;
	}
	
	// calculate moment of inertia coefficient.
	this.iSum = 1/(1/a.i + 1/b.i);
	
	// calculate bias velocity
	var maxBias = this.maxBias;
	this.bias = clamp(-bias_coef(this.errorBias, dt)*pdist/dt, -maxBias, maxBias);
	
	// compute max impulse
	this.jMax = this.maxForce * dt;

	// If the bias is 0, the joint is not at a limit. Reset the impulse.
	if(!this.bias) this.jAcc = 0;
};

RotaryLimitJoint.prototype.applyCachedImpulse = function(dt_coef)
{
	var a = this.a;
	var b = this.b;
	
	var j = this.jAcc*dt_coef;
	a.w -= j*a.i_inv;
	b.w += j*b.i_inv;
};

RotaryLimitJoint.prototype.applyImpulse = function()
{
	if(!this.bias) return; // early exit

	var a = this.a;
	var b = this.b;
	
	// compute relative rotational velocity
	var wr = b.w - a.w;
	
	// compute normal impulse	
	var j = -(this.bias + wr)*this.iSum;
	var jOld = this.jAcc;
	if(this.bias < 0){
		this.jAcc = clamp(jOld + j, 0, this.jMax);
	} else {
		this.jAcc = clamp(jOld + j, -this.jMax, 0);
	}
	j = this.jAcc - jOld;
	
	// apply impulse
	a.w -= j*a.i_inv;
	b.w += j*b.i_inv;
};

RotaryLimitJoint.prototype.getImpulse = function()
{
	return Math.abs(joint.jAcc);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var RatchetJoint = cp.RatchetJoint = function(a, b, phase, ratchet)
{
	Constraint.call(this, a, b);

	this.angle = 0;
	this.phase = phase;
	this.ratchet = ratchet;
	
	// STATIC_BODY_CHECK
	this.angle = (b ? b.a : 0) - (a ? a.a : 0);
	
	this.iSum = this.bias = this.jAcc = this.jMax = 0;
};

RatchetJoint.prototype = Object.create(Constraint.prototype);

RatchetJoint.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	var angle = this.angle;
	var phase = this.phase;
	var ratchet = this.ratchet;
	
	var delta = b.a - a.a;
	var diff = angle - delta;
	var pdist = 0;
	
	if(diff*ratchet > 0){
		pdist = diff;
	} else {
		this.angle = Math.floor((delta - phase)/ratchet)*ratchet + phase;
	}
	
	// calculate moment of inertia coefficient.
	this.iSum = 1/(a.i_inv + b.i_inv);
	
	// calculate bias velocity
	var maxBias = this.maxBias;
	this.bias = clamp(-bias_coef(this.errorBias, dt)*pdist/dt, -maxBias, maxBias);
	
	// compute max impulse
	this.jMax = this.maxForce * dt;

	// If the bias is 0, the joint is not at a limit. Reset the impulse.
	if(!this.bias) this.jAcc = 0;
};

RatchetJoint.prototype.applyCachedImpulse = function(dt_coef)
{
	var a = this.a;
	var b = this.b;
	
	var j = this.jAcc*dt_coef;
	a.w -= j*a.i_inv;
	b.w += j*b.i_inv;
};

RatchetJoint.prototype.applyImpulse = function()
{
	if(!this.bias) return; // early exit

	var a = this.a;
	var b = this.b;
	
	// compute relative rotational velocity
	var wr = b.w - a.w;
	var ratchet = this.ratchet;
	
	// compute normal impulse	
	var j = -(this.bias + wr)*this.iSum;
	var jOld = this.jAcc;
	this.jAcc = clamp((jOld + j)*ratchet, 0, this.jMax*Math.abs(ratchet))/ratchet;
	j = this.jAcc - jOld;
	
	// apply impulse
	a.w -= j*a.i_inv;
	b.w += j*b.i_inv;
};

RatchetJoint.prototype.getImpulse = function(joint)
{
	return Math.abs(joint.jAcc);
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var GearJoint = cp.GearJoint = function(a, b, phase, ratio)
{
	Constraint.call(this, a, b);
	
	this.phase = phase;
	this.ratio = ratio;
	this.ratio_inv = 1/ratio;
	
	this.jAcc = 0;

	this.iSum = this.bias = this.jMax = 0;
};

GearJoint.prototype = Object.create(Constraint.prototype);

GearJoint.prototype.preStep = function(dt)
{
	var a = this.a;
	var b = this.b;
	
	// calculate moment of inertia coefficient.
	this.iSum = 1/(a.i_inv*this.ratio_inv + this.ratio*b.i_inv);
	
	// calculate bias velocity
	var maxBias = this.maxBias;
	this.bias = clamp(-bias_coef(this.errorBias, dt)*(b.a*this.ratio - a.a - this.phase)/dt, -maxBias, maxBias);
	
	// compute max impulse
	this.jMax = this.maxForce * dt;
};

GearJoint.prototype.applyCachedImpulse = function(dt_coef)
{
	var a = this.a;
	var b = this.b;
	
	var j = this.jAcc*dt_coef;
	a.w -= j*a.i_inv*this.ratio_inv;
	b.w += j*b.i_inv;
};

GearJoint.prototype.applyImpulse = function()
{
	var a = this.a;
	var b = this.b;
	
	// compute relative rotational velocity
	var wr = b.w*this.ratio - a.w;
	
	// compute normal impulse	
	var j = (this.bias - wr)*this.iSum;
	var jOld = this.jAcc;
	this.jAcc = clamp(jOld + j, -this.jMax, this.jMax);
	j = this.jAcc - jOld;
	
	// apply impulse
	a.w -= j*a.i_inv*this.ratio_inv;
	b.w += j*b.i_inv;
};

GearJoint.prototype.getImpulse = function()
{
	return Math.abs(this.jAcc);
};

GearJoint.prototype.setRatio = function(value)
{
	this.ratio = value;
	this.ratio_inv = 1/value;
	this.activateBodies();
};

/* Copyright (c) 2007 Scott Lembcke
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

var SimpleMotor = cp.SimpleMotor = function(a, b, rate)
{
	Constraint.call(this, a, b);
	
	this.rate = rate;
	
	this.jAcc = 0;

	this.iSum = this.jMax = 0;
};

SimpleMotor.prototype = Object.create(Constraint.prototype);

SimpleMotor.prototype.preStep = function(dt)
{
	// calculate moment of inertia coefficient.
	this.iSum = 1/(this.a.i_inv + this.b.i_inv);
	
	// compute max impulse
	this.jMax = this.maxForce * dt;
};

SimpleMotor.prototype.applyCachedImpulse = function(dt_coef)
{
	var a = this.a;
	var b = this.b;
	
	var j = this.jAcc*dt_coef;
	a.w -= j*a.i_inv;
	b.w += j*b.i_inv;
};

SimpleMotor.prototype.applyImpulse = function()
{
	var a = this.a;
	var b = this.b;
	
	// compute relative rotational velocity
	var wr = b.w - a.w + this.rate;
	
	// compute normal impulse	
	var j = -wr*this.iSum;
	var jOld = this.jAcc;
	this.jAcc = clamp(jOld + j, -this.jMax, this.jMax);
	j = this.jAcc - jOld;
	
	// apply impulse
	a.w -= j*a.i_inv;
	b.w += j*b.i_inv;
};

SimpleMotor.prototype.getImpulse = function()
{
	return Math.abs(this.jAcc);
};

})();

},{}],2:[function(_dereq_,module,exports){
/**
 * The buffer module from node.js, for the browser.
 *
 * Author:   Feross Aboukhadijeh <feross@feross.org> <http://feross.org>
 * License:  MIT
 *
 * `npm install buffer`
 */

var base64 = _dereq_('base64-js')
var ieee754 = _dereq_('ieee754')

exports.Buffer = Buffer
exports.SlowBuffer = Buffer
exports.INSPECT_MAX_BYTES = 50
Buffer.poolSize = 8192

/**
 * If `Buffer._useTypedArrays`:
 *   === true    Use Uint8Array implementation (fastest)
 *   === false   Use Object implementation (compatible down to IE6)
 */
Buffer._useTypedArrays = (function () {
   // Detect if browser supports Typed Arrays. Supported browsers are IE 10+,
   // Firefox 4+, Chrome 7+, Safari 5.1+, Opera 11.6+, iOS 4.2+.
  if (typeof Uint8Array === 'undefined' || typeof ArrayBuffer === 'undefined')
    return false

  // Does the browser support adding properties to `Uint8Array` instances? If
  // not, then that's the same as no `Uint8Array` support. We need to be able to
  // add all the node Buffer API methods.
  // Relevant Firefox bug: https://bugzilla.mozilla.org/show_bug.cgi?id=695438
  try {
    var arr = new Uint8Array(0)
    arr.foo = function () { return 42 }
    return 42 === arr.foo() &&
        typeof arr.subarray === 'function' // Chrome 9-10 lack `subarray`
  } catch (e) {
    return false
  }
})()

/**
 * Class: Buffer
 * =============
 *
 * The Buffer constructor returns instances of `Uint8Array` that are augmented
 * with function properties for all the node `Buffer` API functions. We use
 * `Uint8Array` so that square bracket notation works as expected -- it returns
 * a single octet.
 *
 * By augmenting the instances, we can avoid modifying the `Uint8Array`
 * prototype.
 */
function Buffer (subject, encoding, noZero) {
  if (!(this instanceof Buffer))
    return new Buffer(subject, encoding, noZero)

  var type = typeof subject

  // Workaround: node's base64 implementation allows for non-padded strings
  // while base64-js does not.
  if (encoding === 'base64' && type === 'string') {
    subject = stringtrim(subject)
    while (subject.length % 4 !== 0) {
      subject = subject + '='
    }
  }

  // Find the length
  var length
  if (type === 'number')
    length = coerce(subject)
  else if (type === 'string')
    length = Buffer.byteLength(subject, encoding)
  else if (type === 'object')
    length = coerce(subject.length) // Assume object is an array
  else
    throw new Error('First argument needs to be a number, array or string.')

  var buf
  if (Buffer._useTypedArrays) {
    // Preferred: Return an augmented `Uint8Array` instance for best performance
    buf = augment(new Uint8Array(length))
  } else {
    // Fallback: Return THIS instance of Buffer (created by `new`)
    buf = this
    buf.length = length
    buf._isBuffer = true
  }

  var i
  if (Buffer._useTypedArrays && typeof Uint8Array === 'function' &&
      subject instanceof Uint8Array) {
    // Speed optimization -- use set if we're copying from a Uint8Array
    buf._set(subject)
  } else if (isArrayish(subject)) {
    // Treat array-ish objects as a byte array
    for (i = 0; i < length; i++) {
      if (Buffer.isBuffer(subject))
        buf[i] = subject.readUInt8(i)
      else
        buf[i] = subject[i]
    }
  } else if (type === 'string') {
    buf.write(subject, 0, encoding)
  } else if (type === 'number' && !Buffer._useTypedArrays && !noZero) {
    for (i = 0; i < length; i++) {
      buf[i] = 0
    }
  }

  return buf
}

// STATIC METHODS
// ==============

Buffer.isEncoding = function (encoding) {
  switch (String(encoding).toLowerCase()) {
    case 'hex':
    case 'utf8':
    case 'utf-8':
    case 'ascii':
    case 'binary':
    case 'base64':
    case 'raw':
    case 'ucs2':
    case 'ucs-2':
    case 'utf16le':
    case 'utf-16le':
      return true
    default:
      return false
  }
}

Buffer.isBuffer = function (b) {
  return !!(b !== null && b !== undefined && b._isBuffer)
}

Buffer.byteLength = function (str, encoding) {
  var ret
  str = str + ''
  switch (encoding || 'utf8') {
    case 'hex':
      ret = str.length / 2
      break
    case 'utf8':
    case 'utf-8':
      ret = utf8ToBytes(str).length
      break
    case 'ascii':
    case 'binary':
    case 'raw':
      ret = str.length
      break
    case 'base64':
      ret = base64ToBytes(str).length
      break
    case 'ucs2':
    case 'ucs-2':
    case 'utf16le':
    case 'utf-16le':
      ret = str.length * 2
      break
    default:
      throw new Error('Unknown encoding')
  }
  return ret
}

Buffer.concat = function (list, totalLength) {
  assert(isArray(list), 'Usage: Buffer.concat(list, [totalLength])\n' +
      'list should be an Array.')

  if (list.length === 0) {
    return new Buffer(0)
  } else if (list.length === 1) {
    return list[0]
  }

  var i
  if (typeof totalLength !== 'number') {
    totalLength = 0
    for (i = 0; i < list.length; i++) {
      totalLength += list[i].length
    }
  }

  var buf = new Buffer(totalLength)
  var pos = 0
  for (i = 0; i < list.length; i++) {
    var item = list[i]
    item.copy(buf, pos)
    pos += item.length
  }
  return buf
}

// BUFFER INSTANCE METHODS
// =======================

function _hexWrite (buf, string, offset, length) {
  offset = Number(offset) || 0
  var remaining = buf.length - offset
  if (!length) {
    length = remaining
  } else {
    length = Number(length)
    if (length > remaining) {
      length = remaining
    }
  }

  // must be an even number of digits
  var strLen = string.length
  assert(strLen % 2 === 0, 'Invalid hex string')

  if (length > strLen / 2) {
    length = strLen / 2
  }
  for (var i = 0; i < length; i++) {
    var byte = parseInt(string.substr(i * 2, 2), 16)
    assert(!isNaN(byte), 'Invalid hex string')
    buf[offset + i] = byte
  }
  Buffer._charsWritten = i * 2
  return i
}

function _utf8Write (buf, string, offset, length) {
  var charsWritten = Buffer._charsWritten =
    blitBuffer(utf8ToBytes(string), buf, offset, length)
  return charsWritten
}

function _asciiWrite (buf, string, offset, length) {
  var charsWritten = Buffer._charsWritten =
    blitBuffer(asciiToBytes(string), buf, offset, length)
  return charsWritten
}

function _binaryWrite (buf, string, offset, length) {
  return _asciiWrite(buf, string, offset, length)
}

function _base64Write (buf, string, offset, length) {
  var charsWritten = Buffer._charsWritten =
    blitBuffer(base64ToBytes(string), buf, offset, length)
  return charsWritten
}

function _utf16leWrite (buf, string, offset, length) {
  var charsWritten = Buffer._charsWritten =
    blitBuffer(utf16leToBytes(string), buf, offset, length)
  return charsWritten
}

Buffer.prototype.write = function (string, offset, length, encoding) {
  // Support both (string, offset, length, encoding)
  // and the legacy (string, encoding, offset, length)
  if (isFinite(offset)) {
    if (!isFinite(length)) {
      encoding = length
      length = undefined
    }
  } else {  // legacy
    var swap = encoding
    encoding = offset
    offset = length
    length = swap
  }

  offset = Number(offset) || 0
  var remaining = this.length - offset
  if (!length) {
    length = remaining
  } else {
    length = Number(length)
    if (length > remaining) {
      length = remaining
    }
  }
  encoding = String(encoding || 'utf8').toLowerCase()

  var ret
  switch (encoding) {
    case 'hex':
      ret = _hexWrite(this, string, offset, length)
      break
    case 'utf8':
    case 'utf-8':
      ret = _utf8Write(this, string, offset, length)
      break
    case 'ascii':
      ret = _asciiWrite(this, string, offset, length)
      break
    case 'binary':
      ret = _binaryWrite(this, string, offset, length)
      break
    case 'base64':
      ret = _base64Write(this, string, offset, length)
      break
    case 'ucs2':
    case 'ucs-2':
    case 'utf16le':
    case 'utf-16le':
      ret = _utf16leWrite(this, string, offset, length)
      break
    default:
      throw new Error('Unknown encoding')
  }
  return ret
}

Buffer.prototype.toString = function (encoding, start, end) {
  var self = this

  encoding = String(encoding || 'utf8').toLowerCase()
  start = Number(start) || 0
  end = (end !== undefined)
    ? Number(end)
    : end = self.length

  // Fastpath empty strings
  if (end === start)
    return ''

  var ret
  switch (encoding) {
    case 'hex':
      ret = _hexSlice(self, start, end)
      break
    case 'utf8':
    case 'utf-8':
      ret = _utf8Slice(self, start, end)
      break
    case 'ascii':
      ret = _asciiSlice(self, start, end)
      break
    case 'binary':
      ret = _binarySlice(self, start, end)
      break
    case 'base64':
      ret = _base64Slice(self, start, end)
      break
    case 'ucs2':
    case 'ucs-2':
    case 'utf16le':
    case 'utf-16le':
      ret = _utf16leSlice(self, start, end)
      break
    default:
      throw new Error('Unknown encoding')
  }
  return ret
}

Buffer.prototype.toJSON = function () {
  return {
    type: 'Buffer',
    data: Array.prototype.slice.call(this._arr || this, 0)
  }
}

// copy(targetBuffer, targetStart=0, sourceStart=0, sourceEnd=buffer.length)
Buffer.prototype.copy = function (target, target_start, start, end) {
  var source = this

  if (!start) start = 0
  if (!end && end !== 0) end = this.length
  if (!target_start) target_start = 0

  // Copy 0 bytes; we're done
  if (end === start) return
  if (target.length === 0 || source.length === 0) return

  // Fatal error conditions
  assert(end >= start, 'sourceEnd < sourceStart')
  assert(target_start >= 0 && target_start < target.length,
      'targetStart out of bounds')
  assert(start >= 0 && start < source.length, 'sourceStart out of bounds')
  assert(end >= 0 && end <= source.length, 'sourceEnd out of bounds')

  // Are we oob?
  if (end > this.length)
    end = this.length
  if (target.length - target_start < end - start)
    end = target.length - target_start + start

  // copy!
  for (var i = 0; i < end - start; i++)
    target[i + target_start] = this[i + start]
}

function _base64Slice (buf, start, end) {
  if (start === 0 && end === buf.length) {
    return base64.fromByteArray(buf)
  } else {
    return base64.fromByteArray(buf.slice(start, end))
  }
}

function _utf8Slice (buf, start, end) {
  var res = ''
  var tmp = ''
  end = Math.min(buf.length, end)

  for (var i = start; i < end; i++) {
    if (buf[i] <= 0x7F) {
      res += decodeUtf8Char(tmp) + String.fromCharCode(buf[i])
      tmp = ''
    } else {
      tmp += '%' + buf[i].toString(16)
    }
  }

  return res + decodeUtf8Char(tmp)
}

function _asciiSlice (buf, start, end) {
  var ret = ''
  end = Math.min(buf.length, end)

  for (var i = start; i < end; i++)
    ret += String.fromCharCode(buf[i])
  return ret
}

function _binarySlice (buf, start, end) {
  return _asciiSlice(buf, start, end)
}

function _hexSlice (buf, start, end) {
  var len = buf.length

  if (!start || start < 0) start = 0
  if (!end || end < 0 || end > len) end = len

  var out = ''
  for (var i = start; i < end; i++) {
    out += toHex(buf[i])
  }
  return out
}

function _utf16leSlice (buf, start, end) {
  var bytes = buf.slice(start, end)
  var res = ''
  for (var i = 0; i < bytes.length; i += 2) {
    res += String.fromCharCode(bytes[i] + bytes[i+1] * 256)
  }
  return res
}

Buffer.prototype.slice = function (start, end) {
  var len = this.length
  start = clamp(start, len, 0)
  end = clamp(end, len, len)

  if (Buffer._useTypedArrays) {
    return augment(this.subarray(start, end))
  } else {
    var sliceLen = end - start
    var newBuf = new Buffer(sliceLen, undefined, true)
    for (var i = 0; i < sliceLen; i++) {
      newBuf[i] = this[i + start]
    }
    return newBuf
  }
}

// `get` will be removed in Node 0.13+
Buffer.prototype.get = function (offset) {
  console.log('.get() is deprecated. Access using array indexes instead.')
  return this.readUInt8(offset)
}

// `set` will be removed in Node 0.13+
Buffer.prototype.set = function (v, offset) {
  console.log('.set() is deprecated. Access using array indexes instead.')
  return this.writeUInt8(v, offset)
}

Buffer.prototype.readUInt8 = function (offset, noAssert) {
  if (!noAssert) {
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset < this.length, 'Trying to read beyond buffer length')
  }

  if (offset >= this.length)
    return

  return this[offset]
}

function _readUInt16 (buf, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 1 < buf.length, 'Trying to read beyond buffer length')
  }

  var len = buf.length
  if (offset >= len)
    return

  var val
  if (littleEndian) {
    val = buf[offset]
    if (offset + 1 < len)
      val |= buf[offset + 1] << 8
  } else {
    val = buf[offset] << 8
    if (offset + 1 < len)
      val |= buf[offset + 1]
  }
  return val
}

Buffer.prototype.readUInt16LE = function (offset, noAssert) {
  return _readUInt16(this, offset, true, noAssert)
}

Buffer.prototype.readUInt16BE = function (offset, noAssert) {
  return _readUInt16(this, offset, false, noAssert)
}

function _readUInt32 (buf, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 3 < buf.length, 'Trying to read beyond buffer length')
  }

  var len = buf.length
  if (offset >= len)
    return

  var val
  if (littleEndian) {
    if (offset + 2 < len)
      val = buf[offset + 2] << 16
    if (offset + 1 < len)
      val |= buf[offset + 1] << 8
    val |= buf[offset]
    if (offset + 3 < len)
      val = val + (buf[offset + 3] << 24 >>> 0)
  } else {
    if (offset + 1 < len)
      val = buf[offset + 1] << 16
    if (offset + 2 < len)
      val |= buf[offset + 2] << 8
    if (offset + 3 < len)
      val |= buf[offset + 3]
    val = val + (buf[offset] << 24 >>> 0)
  }
  return val
}

Buffer.prototype.readUInt32LE = function (offset, noAssert) {
  return _readUInt32(this, offset, true, noAssert)
}

Buffer.prototype.readUInt32BE = function (offset, noAssert) {
  return _readUInt32(this, offset, false, noAssert)
}

Buffer.prototype.readInt8 = function (offset, noAssert) {
  if (!noAssert) {
    assert(offset !== undefined && offset !== null,
        'missing offset')
    assert(offset < this.length, 'Trying to read beyond buffer length')
  }

  if (offset >= this.length)
    return

  var neg = this[offset] & 0x80
  if (neg)
    return (0xff - this[offset] + 1) * -1
  else
    return this[offset]
}

function _readInt16 (buf, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 1 < buf.length, 'Trying to read beyond buffer length')
  }

  var len = buf.length
  if (offset >= len)
    return

  var val = _readUInt16(buf, offset, littleEndian, true)
  var neg = val & 0x8000
  if (neg)
    return (0xffff - val + 1) * -1
  else
    return val
}

Buffer.prototype.readInt16LE = function (offset, noAssert) {
  return _readInt16(this, offset, true, noAssert)
}

Buffer.prototype.readInt16BE = function (offset, noAssert) {
  return _readInt16(this, offset, false, noAssert)
}

function _readInt32 (buf, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 3 < buf.length, 'Trying to read beyond buffer length')
  }

  var len = buf.length
  if (offset >= len)
    return

  var val = _readUInt32(buf, offset, littleEndian, true)
  var neg = val & 0x80000000
  if (neg)
    return (0xffffffff - val + 1) * -1
  else
    return val
}

Buffer.prototype.readInt32LE = function (offset, noAssert) {
  return _readInt32(this, offset, true, noAssert)
}

Buffer.prototype.readInt32BE = function (offset, noAssert) {
  return _readInt32(this, offset, false, noAssert)
}

function _readFloat (buf, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset + 3 < buf.length, 'Trying to read beyond buffer length')
  }

  return ieee754.read(buf, offset, littleEndian, 23, 4)
}

Buffer.prototype.readFloatLE = function (offset, noAssert) {
  return _readFloat(this, offset, true, noAssert)
}

Buffer.prototype.readFloatBE = function (offset, noAssert) {
  return _readFloat(this, offset, false, noAssert)
}

function _readDouble (buf, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset + 7 < buf.length, 'Trying to read beyond buffer length')
  }

  return ieee754.read(buf, offset, littleEndian, 52, 8)
}

Buffer.prototype.readDoubleLE = function (offset, noAssert) {
  return _readDouble(this, offset, true, noAssert)
}

Buffer.prototype.readDoubleBE = function (offset, noAssert) {
  return _readDouble(this, offset, false, noAssert)
}

Buffer.prototype.writeUInt8 = function (value, offset, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset < this.length, 'trying to write beyond buffer length')
    verifuint(value, 0xff)
  }

  if (offset >= this.length) return

  this[offset] = value
}

function _writeUInt16 (buf, value, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 1 < buf.length, 'trying to write beyond buffer length')
    verifuint(value, 0xffff)
  }

  var len = buf.length
  if (offset >= len)
    return

  for (var i = 0, j = Math.min(len - offset, 2); i < j; i++) {
    buf[offset + i] =
        (value & (0xff << (8 * (littleEndian ? i : 1 - i)))) >>>
            (littleEndian ? i : 1 - i) * 8
  }
}

Buffer.prototype.writeUInt16LE = function (value, offset, noAssert) {
  _writeUInt16(this, value, offset, true, noAssert)
}

Buffer.prototype.writeUInt16BE = function (value, offset, noAssert) {
  _writeUInt16(this, value, offset, false, noAssert)
}

function _writeUInt32 (buf, value, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 3 < buf.length, 'trying to write beyond buffer length')
    verifuint(value, 0xffffffff)
  }

  var len = buf.length
  if (offset >= len)
    return

  for (var i = 0, j = Math.min(len - offset, 4); i < j; i++) {
    buf[offset + i] =
        (value >>> (littleEndian ? i : 3 - i) * 8) & 0xff
  }
}

Buffer.prototype.writeUInt32LE = function (value, offset, noAssert) {
  _writeUInt32(this, value, offset, true, noAssert)
}

Buffer.prototype.writeUInt32BE = function (value, offset, noAssert) {
  _writeUInt32(this, value, offset, false, noAssert)
}

Buffer.prototype.writeInt8 = function (value, offset, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset < this.length, 'Trying to write beyond buffer length')
    verifsint(value, 0x7f, -0x80)
  }

  if (offset >= this.length)
    return

  if (value >= 0)
    this.writeUInt8(value, offset, noAssert)
  else
    this.writeUInt8(0xff + value + 1, offset, noAssert)
}

function _writeInt16 (buf, value, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 1 < buf.length, 'Trying to write beyond buffer length')
    verifsint(value, 0x7fff, -0x8000)
  }

  var len = buf.length
  if (offset >= len)
    return

  if (value >= 0)
    _writeUInt16(buf, value, offset, littleEndian, noAssert)
  else
    _writeUInt16(buf, 0xffff + value + 1, offset, littleEndian, noAssert)
}

Buffer.prototype.writeInt16LE = function (value, offset, noAssert) {
  _writeInt16(this, value, offset, true, noAssert)
}

Buffer.prototype.writeInt16BE = function (value, offset, noAssert) {
  _writeInt16(this, value, offset, false, noAssert)
}

function _writeInt32 (buf, value, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 3 < buf.length, 'Trying to write beyond buffer length')
    verifsint(value, 0x7fffffff, -0x80000000)
  }

  var len = buf.length
  if (offset >= len)
    return

  if (value >= 0)
    _writeUInt32(buf, value, offset, littleEndian, noAssert)
  else
    _writeUInt32(buf, 0xffffffff + value + 1, offset, littleEndian, noAssert)
}

Buffer.prototype.writeInt32LE = function (value, offset, noAssert) {
  _writeInt32(this, value, offset, true, noAssert)
}

Buffer.prototype.writeInt32BE = function (value, offset, noAssert) {
  _writeInt32(this, value, offset, false, noAssert)
}

function _writeFloat (buf, value, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 3 < buf.length, 'Trying to write beyond buffer length')
    verifIEEE754(value, 3.4028234663852886e+38, -3.4028234663852886e+38)
  }

  var len = buf.length
  if (offset >= len)
    return

  ieee754.write(buf, value, offset, littleEndian, 23, 4)
}

Buffer.prototype.writeFloatLE = function (value, offset, noAssert) {
  _writeFloat(this, value, offset, true, noAssert)
}

Buffer.prototype.writeFloatBE = function (value, offset, noAssert) {
  _writeFloat(this, value, offset, false, noAssert)
}

function _writeDouble (buf, value, offset, littleEndian, noAssert) {
  if (!noAssert) {
    assert(value !== undefined && value !== null, 'missing value')
    assert(typeof littleEndian === 'boolean', 'missing or invalid endian')
    assert(offset !== undefined && offset !== null, 'missing offset')
    assert(offset + 7 < buf.length,
        'Trying to write beyond buffer length')
    verifIEEE754(value, 1.7976931348623157E+308, -1.7976931348623157E+308)
  }

  var len = buf.length
  if (offset >= len)
    return

  ieee754.write(buf, value, offset, littleEndian, 52, 8)
}

Buffer.prototype.writeDoubleLE = function (value, offset, noAssert) {
  _writeDouble(this, value, offset, true, noAssert)
}

Buffer.prototype.writeDoubleBE = function (value, offset, noAssert) {
  _writeDouble(this, value, offset, false, noAssert)
}

// fill(value, start=0, end=buffer.length)
Buffer.prototype.fill = function (value, start, end) {
  if (!value) value = 0
  if (!start) start = 0
  if (!end) end = this.length

  if (typeof value === 'string') {
    value = value.charCodeAt(0)
  }

  assert(typeof value === 'number' && !isNaN(value), 'value is not a number')
  assert(end >= start, 'end < start')

  // Fill 0 bytes; we're done
  if (end === start) return
  if (this.length === 0) return

  assert(start >= 0 && start < this.length, 'start out of bounds')
  assert(end >= 0 && end <= this.length, 'end out of bounds')

  for (var i = start; i < end; i++) {
    this[i] = value
  }
}

Buffer.prototype.inspect = function () {
  var out = []
  var len = this.length
  for (var i = 0; i < len; i++) {
    out[i] = toHex(this[i])
    if (i === exports.INSPECT_MAX_BYTES) {
      out[i + 1] = '...'
      break
    }
  }
  return '<Buffer ' + out.join(' ') + '>'
}

/**
 * Creates a new `ArrayBuffer` with the *copied* memory of the buffer instance.
 * Added in Node 0.12. Only available in browsers that support ArrayBuffer.
 */
Buffer.prototype.toArrayBuffer = function () {
  if (typeof Uint8Array === 'function') {
    if (Buffer._useTypedArrays) {
      return (new Buffer(this)).buffer
    } else {
      var buf = new Uint8Array(this.length)
      for (var i = 0, len = buf.length; i < len; i += 1)
        buf[i] = this[i]
      return buf.buffer
    }
  } else {
    throw new Error('Buffer.toArrayBuffer not supported in this browser')
  }
}

// HELPER FUNCTIONS
// ================

function stringtrim (str) {
  if (str.trim) return str.trim()
  return str.replace(/^\s+|\s+$/g, '')
}

var BP = Buffer.prototype

/**
 * Augment the Uint8Array *instance* (not the class!) with Buffer methods
 */
function augment (arr) {
  arr._isBuffer = true

  // save reference to original Uint8Array get/set methods before overwriting
  arr._get = arr.get
  arr._set = arr.set

  // deprecated, will be removed in node 0.13+
  arr.get = BP.get
  arr.set = BP.set

  arr.write = BP.write
  arr.toString = BP.toString
  arr.toLocaleString = BP.toString
  arr.toJSON = BP.toJSON
  arr.copy = BP.copy
  arr.slice = BP.slice
  arr.readUInt8 = BP.readUInt8
  arr.readUInt16LE = BP.readUInt16LE
  arr.readUInt16BE = BP.readUInt16BE
  arr.readUInt32LE = BP.readUInt32LE
  arr.readUInt32BE = BP.readUInt32BE
  arr.readInt8 = BP.readInt8
  arr.readInt16LE = BP.readInt16LE
  arr.readInt16BE = BP.readInt16BE
  arr.readInt32LE = BP.readInt32LE
  arr.readInt32BE = BP.readInt32BE
  arr.readFloatLE = BP.readFloatLE
  arr.readFloatBE = BP.readFloatBE
  arr.readDoubleLE = BP.readDoubleLE
  arr.readDoubleBE = BP.readDoubleBE
  arr.writeUInt8 = BP.writeUInt8
  arr.writeUInt16LE = BP.writeUInt16LE
  arr.writeUInt16BE = BP.writeUInt16BE
  arr.writeUInt32LE = BP.writeUInt32LE
  arr.writeUInt32BE = BP.writeUInt32BE
  arr.writeInt8 = BP.writeInt8
  arr.writeInt16LE = BP.writeInt16LE
  arr.writeInt16BE = BP.writeInt16BE
  arr.writeInt32LE = BP.writeInt32LE
  arr.writeInt32BE = BP.writeInt32BE
  arr.writeFloatLE = BP.writeFloatLE
  arr.writeFloatBE = BP.writeFloatBE
  arr.writeDoubleLE = BP.writeDoubleLE
  arr.writeDoubleBE = BP.writeDoubleBE
  arr.fill = BP.fill
  arr.inspect = BP.inspect
  arr.toArrayBuffer = BP.toArrayBuffer

  return arr
}

// slice(start, end)
function clamp (index, len, defaultValue) {
  if (typeof index !== 'number') return defaultValue
  index = ~~index;  // Coerce to integer.
  if (index >= len) return len
  if (index >= 0) return index
  index += len
  if (index >= 0) return index
  return 0
}

function coerce (length) {
  // Coerce length to a number (possibly NaN), round up
  // in case it's fractional (e.g. 123.456) then do a
  // double negate to coerce a NaN to 0. Easy, right?
  length = ~~Math.ceil(+length)
  return length < 0 ? 0 : length
}

function isArray (subject) {
  return (Array.isArray || function (subject) {
    return Object.prototype.toString.call(subject) === '[object Array]'
  })(subject)
}

function isArrayish (subject) {
  return isArray(subject) || Buffer.isBuffer(subject) ||
      subject && typeof subject === 'object' &&
      typeof subject.length === 'number'
}

function toHex (n) {
  if (n < 16) return '0' + n.toString(16)
  return n.toString(16)
}

function utf8ToBytes (str) {
  var byteArray = []
  for (var i = 0; i < str.length; i++) {
    var b = str.charCodeAt(i)
    if (b <= 0x7F)
      byteArray.push(str.charCodeAt(i))
    else {
      var start = i
      if (b >= 0xD800 && b <= 0xDFFF) i++
      var h = encodeURIComponent(str.slice(start, i+1)).substr(1).split('%')
      for (var j = 0; j < h.length; j++)
        byteArray.push(parseInt(h[j], 16))
    }
  }
  return byteArray
}

function asciiToBytes (str) {
  var byteArray = []
  for (var i = 0; i < str.length; i++) {
    // Node's code seems to be doing this and not & 0x7F..
    byteArray.push(str.charCodeAt(i) & 0xFF)
  }
  return byteArray
}

function utf16leToBytes (str) {
  var c, hi, lo
  var byteArray = []
  for (var i = 0; i < str.length; i++) {
    c = str.charCodeAt(i)
    hi = c >> 8
    lo = c % 256
    byteArray.push(lo)
    byteArray.push(hi)
  }

  return byteArray
}

function base64ToBytes (str) {
  return base64.toByteArray(str)
}

function blitBuffer (src, dst, offset, length) {
  var pos
  for (var i = 0; i < length; i++) {
    if ((i + offset >= dst.length) || (i >= src.length))
      break
    dst[i + offset] = src[i]
  }
  return i
}

function decodeUtf8Char (str) {
  try {
    return decodeURIComponent(str)
  } catch (err) {
    return String.fromCharCode(0xFFFD) // UTF 8 invalid char
  }
}

/*
 * We have to make sure that the value is a valid integer. This means that it
 * is non-negative. It has no fractional component and that it does not
 * exceed the maximum allowed value.
 */
function verifuint (value, max) {
  assert(typeof value === 'number', 'cannot write a non-number as a number')
  assert(value >= 0,
      'specified a negative value for writing an unsigned value')
  assert(value <= max, 'value is larger than maximum value for type')
  assert(Math.floor(value) === value, 'value has a fractional component')
}

function verifsint (value, max, min) {
  assert(typeof value === 'number', 'cannot write a non-number as a number')
  assert(value <= max, 'value larger than maximum allowed value')
  assert(value >= min, 'value smaller than minimum allowed value')
  assert(Math.floor(value) === value, 'value has a fractional component')
}

function verifIEEE754 (value, max, min) {
  assert(typeof value === 'number', 'cannot write a non-number as a number')
  assert(value <= max, 'value larger than maximum allowed value')
  assert(value >= min, 'value smaller than minimum allowed value')
}

function assert (test, message) {
  if (!test) throw new Error(message || 'Failed assertion')
}

},{"base64-js":3,"ieee754":4}],3:[function(_dereq_,module,exports){
var lookup = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';

;(function (exports) {
	'use strict';

  var Arr = (typeof Uint8Array !== 'undefined')
    ? Uint8Array
    : Array

	var ZERO   = '0'.charCodeAt(0)
	var PLUS   = '+'.charCodeAt(0)
	var SLASH  = '/'.charCodeAt(0)
	var NUMBER = '0'.charCodeAt(0)
	var LOWER  = 'a'.charCodeAt(0)
	var UPPER  = 'A'.charCodeAt(0)

	function decode (elt) {
		var code = elt.charCodeAt(0)
		if (code === PLUS)
			return 62 // '+'
		if (code === SLASH)
			return 63 // '/'
		if (code < NUMBER)
			return -1 //no match
		if (code < NUMBER + 10)
			return code - NUMBER + 26 + 26
		if (code < UPPER + 26)
			return code - UPPER
		if (code < LOWER + 26)
			return code - LOWER + 26
	}

	function b64ToByteArray (b64) {
		var i, j, l, tmp, placeHolders, arr

		if (b64.length % 4 > 0) {
			throw new Error('Invalid string. Length must be a multiple of 4')
		}

		// the number of equal signs (place holders)
		// if there are two placeholders, than the two characters before it
		// represent one byte
		// if there is only one, then the three characters before it represent 2 bytes
		// this is just a cheap hack to not do indexOf twice
		var len = b64.length
		placeHolders = '=' === b64.charAt(len - 2) ? 2 : '=' === b64.charAt(len - 1) ? 1 : 0

		// base64 is 4/3 + up to two characters of the original data
		arr = new Arr(b64.length * 3 / 4 - placeHolders)

		// if there are placeholders, only get up to the last complete 4 chars
		l = placeHolders > 0 ? b64.length - 4 : b64.length

		var L = 0

		function push (v) {
			arr[L++] = v
		}

		for (i = 0, j = 0; i < l; i += 4, j += 3) {
			tmp = (decode(b64.charAt(i)) << 18) | (decode(b64.charAt(i + 1)) << 12) | (decode(b64.charAt(i + 2)) << 6) | decode(b64.charAt(i + 3))
			push((tmp & 0xFF0000) >> 16)
			push((tmp & 0xFF00) >> 8)
			push(tmp & 0xFF)
		}

		if (placeHolders === 2) {
			tmp = (decode(b64.charAt(i)) << 2) | (decode(b64.charAt(i + 1)) >> 4)
			push(tmp & 0xFF)
		} else if (placeHolders === 1) {
			tmp = (decode(b64.charAt(i)) << 10) | (decode(b64.charAt(i + 1)) << 4) | (decode(b64.charAt(i + 2)) >> 2)
			push((tmp >> 8) & 0xFF)
			push(tmp & 0xFF)
		}

		return arr
	}

	function uint8ToBase64 (uint8) {
		var i,
			extraBytes = uint8.length % 3, // if we have 1 byte left, pad 2 bytes
			output = "",
			temp, length

		function encode (num) {
			return lookup.charAt(num)
		}

		function tripletToBase64 (num) {
			return encode(num >> 18 & 0x3F) + encode(num >> 12 & 0x3F) + encode(num >> 6 & 0x3F) + encode(num & 0x3F)
		}

		// go through the array every three bytes, we'll deal with trailing stuff later
		for (i = 0, length = uint8.length - extraBytes; i < length; i += 3) {
			temp = (uint8[i] << 16) + (uint8[i + 1] << 8) + (uint8[i + 2])
			output += tripletToBase64(temp)
		}

		// pad the end with zeros, but make sure to not forget the extra bytes
		switch (extraBytes) {
			case 1:
				temp = uint8[uint8.length - 1]
				output += encode(temp >> 2)
				output += encode((temp << 4) & 0x3F)
				output += '=='
				break
			case 2:
				temp = (uint8[uint8.length - 2] << 8) + (uint8[uint8.length - 1])
				output += encode(temp >> 10)
				output += encode((temp >> 4) & 0x3F)
				output += encode((temp << 2) & 0x3F)
				output += '='
				break
		}

		return output
	}

	module.exports.toByteArray = b64ToByteArray
	module.exports.fromByteArray = uint8ToBase64
}())

},{}],4:[function(_dereq_,module,exports){
exports.read = function(buffer, offset, isLE, mLen, nBytes) {
  var e, m,
      eLen = nBytes * 8 - mLen - 1,
      eMax = (1 << eLen) - 1,
      eBias = eMax >> 1,
      nBits = -7,
      i = isLE ? (nBytes - 1) : 0,
      d = isLE ? -1 : 1,
      s = buffer[offset + i];

  i += d;

  e = s & ((1 << (-nBits)) - 1);
  s >>= (-nBits);
  nBits += eLen;
  for (; nBits > 0; e = e * 256 + buffer[offset + i], i += d, nBits -= 8);

  m = e & ((1 << (-nBits)) - 1);
  e >>= (-nBits);
  nBits += mLen;
  for (; nBits > 0; m = m * 256 + buffer[offset + i], i += d, nBits -= 8);

  if (e === 0) {
    e = 1 - eBias;
  } else if (e === eMax) {
    return m ? NaN : ((s ? -1 : 1) * Infinity);
  } else {
    m = m + Math.pow(2, mLen);
    e = e - eBias;
  }
  return (s ? -1 : 1) * m * Math.pow(2, e - mLen);
};

exports.write = function(buffer, value, offset, isLE, mLen, nBytes) {
  var e, m, c,
      eLen = nBytes * 8 - mLen - 1,
      eMax = (1 << eLen) - 1,
      eBias = eMax >> 1,
      rt = (mLen === 23 ? Math.pow(2, -24) - Math.pow(2, -77) : 0),
      i = isLE ? 0 : (nBytes - 1),
      d = isLE ? 1 : -1,
      s = value < 0 || (value === 0 && 1 / value < 0) ? 1 : 0;

  value = Math.abs(value);

  if (isNaN(value) || value === Infinity) {
    m = isNaN(value) ? 1 : 0;
    e = eMax;
  } else {
    e = Math.floor(Math.log(value) / Math.LN2);
    if (value * (c = Math.pow(2, -e)) < 1) {
      e--;
      c *= 2;
    }
    if (e + eBias >= 1) {
      value += rt / c;
    } else {
      value += rt * Math.pow(2, 1 - eBias);
    }
    if (value * c >= 2) {
      e++;
      c /= 2;
    }

    if (e + eBias >= eMax) {
      m = 0;
      e = eMax;
    } else if (e + eBias >= 1) {
      m = (value * c - 1) * Math.pow(2, mLen);
      e = e + eBias;
    } else {
      m = value * Math.pow(2, eBias - 1) * Math.pow(2, mLen);
      e = 0;
    }
  }

  for (; mLen >= 8; buffer[offset + i] = m & 0xff, i += d, m /= 256, mLen -= 8);

  e = (e << mLen) | m;
  eLen += mLen;
  for (; eLen > 0; buffer[offset + i] = e & 0xff, i += d, e /= 256, eLen -= 8);

  buffer[offset + i - d] |= s * 128;
};

},{}],5:[function(_dereq_,module,exports){
// shim for using process in browser

var process = module.exports = {};

process.nextTick = (function () {
    var canSetImmediate = typeof window !== 'undefined'
    && window.setImmediate;
    var canPost = typeof window !== 'undefined'
    && window.postMessage && window.addEventListener
    ;

    if (canSetImmediate) {
        return function (f) { return window.setImmediate(f) };
    }

    if (canPost) {
        var queue = [];
        window.addEventListener('message', function (ev) {
            var source = ev.source;
            if ((source === window || source === null) && ev.data === 'process-tick') {
                ev.stopPropagation();
                if (queue.length > 0) {
                    var fn = queue.shift();
                    fn();
                }
            }
        }, true);

        return function nextTick(fn) {
            queue.push(fn);
            window.postMessage('process-tick', '*');
        };
    }

    return function nextTick(fn) {
        setTimeout(fn, 0);
    };
})();

process.title = 'browser';
process.browser = true;
process.env = {};
process.argv = [];

process.binding = function (name) {
    throw new Error('process.binding is not supported');
}

// TODO(shtylman)
process.cwd = function () { return '/' };
process.chdir = function (dir) {
    throw new Error('process.chdir is not supported');
};

},{}],6:[function(_dereq_,module,exports){
/**
 * @license
 * pixi.js - v1.5.1
 * Copyright (c) 2012-2014, Mat Groves
 * http://goodboydigital.com/
 *
 * Compiled: 2014-02-13
 *
 * pixi.js is licensed under the MIT License.
 * http://www.opensource.org/licenses/mit-license.php
 */
(function(){var c=this,d=d||{};d.WEBGL_RENDERER=0,d.CANVAS_RENDERER=1,d.VERSION="v1.5.1",d.blendModes={NORMAL:0,ADD:1,MULTIPLY:2,SCREEN:3,OVERLAY:4,DARKEN:5,LIGHTEN:6,COLOR_DODGE:7,COLOR_BURN:8,HARD_LIGHT:9,SOFT_LIGHT:10,DIFFERENCE:11,EXCLUSION:12,HUE:13,SATURATION:14,COLOR:15,LUMINOSITY:16},d.scaleModes={DEFAULT:0,LINEAR:0,NEAREST:1},d.INTERACTION_FREQUENCY=30,d.AUTO_PREVENT_DEFAULT=!0,d.RAD_TO_DEG=180/Math.PI,d.DEG_TO_RAD=Math.PI/180,d.Point=function(a,b){this.x=a||0,this.y=b||0},d.Point.prototype.clone=function(){return new d.Point(this.x,this.y)},d.Point.prototype.constructor=d.Point,d.Point.prototype.set=function(a,b){this.x=a||0,this.y=b||(0!==b?this.x:0)},d.Rectangle=function(a,b,c,d){this.x=a||0,this.y=b||0,this.width=c||0,this.height=d||0},d.Rectangle.prototype.clone=function(){return new d.Rectangle(this.x,this.y,this.width,this.height)},d.Rectangle.prototype.contains=function(a,b){if(this.width<=0||this.height<=0)return!1;var c=this.x;if(a>=c&&a<=c+this.width){var d=this.y;if(b>=d&&b<=d+this.height)return!0}return!1},d.Rectangle.prototype.constructor=d.Rectangle,d.EmptyRectangle=new d.Rectangle(0,0,0,0),d.Polygon=function(a){if(a instanceof Array||(a=Array.prototype.slice.call(arguments)),"number"==typeof a[0]){for(var b=[],c=0,e=a.length;e>c;c+=2)b.push(new d.Point(a[c],a[c+1]));a=b}this.points=a},d.Polygon.prototype.clone=function(){for(var a=[],b=0;b<this.points.length;b++)a.push(this.points[b].clone());return new d.Polygon(a)},d.Polygon.prototype.contains=function(a,b){for(var c=!1,d=0,e=this.points.length-1;d<this.points.length;e=d++){var f=this.points[d].x,g=this.points[d].y,h=this.points[e].x,i=this.points[e].y,j=g>b!=i>b&&(h-f)*(b-g)/(i-g)+f>a;j&&(c=!c)}return c},d.Polygon.prototype.constructor=d.Polygon,d.Circle=function(a,b,c){this.x=a||0,this.y=b||0,this.radius=c||0},d.Circle.prototype.clone=function(){return new d.Circle(this.x,this.y,this.radius)},d.Circle.prototype.contains=function(a,b){if(this.radius<=0)return!1;var c=this.x-a,d=this.y-b,e=this.radius*this.radius;return c*=c,d*=d,e>=c+d},d.Circle.prototype.constructor=d.Circle,d.Ellipse=function(a,b,c,d){this.x=a||0,this.y=b||0,this.width=c||0,this.height=d||0},d.Ellipse.prototype.clone=function(){return new d.Ellipse(this.x,this.y,this.width,this.height)},d.Ellipse.prototype.contains=function(a,b){if(this.width<=0||this.height<=0)return!1;var c=(a-this.x)/this.width,d=(b-this.y)/this.height;return c*=c,d*=d,1>=c+d},d.Ellipse.prototype.getBounds=function(){return new d.Rectangle(this.x,this.y,this.width,this.height)},d.Ellipse.prototype.constructor=d.Ellipse,d.determineMatrixArrayType=function(){return"undefined"!=typeof Float32Array?Float32Array:Array},d.Matrix2=d.determineMatrixArrayType(),d.Matrix=function(){this.a=1,this.b=0,this.c=0,this.d=1,this.tx=0,this.ty=0},d.Matrix.prototype.fromArray=function(a){this.a=a[0],this.b=a[1],this.c=a[3],this.d=a[4],this.tx=a[2],this.ty=a[5]},d.Matrix.prototype.toArray=function(a){this.array||(this.array=new Float32Array(9));var b=this.array;return a?(this.array[0]=this.a,this.array[1]=this.c,this.array[2]=0,this.array[3]=this.b,this.array[4]=this.d,this.array[5]=0,this.array[6]=this.tx,this.array[7]=this.ty,this.array[8]=1):(this.array[0]=this.a,this.array[1]=this.b,this.array[2]=this.tx,this.array[3]=this.c,this.array[4]=this.d,this.array[5]=this.ty,this.array[6]=0,this.array[7]=0,this.array[8]=1),b},d.identityMatrix=new d.Matrix,d.DisplayObject=function(){this.position=new d.Point,this.scale=new d.Point(1,1),this.pivot=new d.Point(0,0),this.rotation=0,this.alpha=1,this.visible=!0,this.hitArea=null,this.buttonMode=!1,this.renderable=!1,this.parent=null,this.stage=null,this.worldAlpha=1,this._interactive=!1,this.defaultCursor="pointer",this.worldTransform=new d.Matrix,this.color=[],this.dynamic=!0,this._sr=0,this._cr=1,this.filterArea=new d.Rectangle(0,0,1,1),this._bounds=new d.Rectangle(0,0,1,1),this._currentBounds=null,this._mask=null},d.DisplayObject.prototype.constructor=d.DisplayObject,d.DisplayObject.prototype.setInteractive=function(a){this.interactive=a},Object.defineProperty(d.DisplayObject.prototype,"interactive",{get:function(){return this._interactive},set:function(a){this._interactive=a,this.stage&&(this.stage.dirty=!0)}}),Object.defineProperty(d.DisplayObject.prototype,"worldVisible",{get:function(){var a=this;do{if(!a.visible)return!1;a=a.parent}while(a);return!0}}),Object.defineProperty(d.DisplayObject.prototype,"mask",{get:function(){return this._mask},set:function(a){this._mask&&(this._mask.isMask=!1),this._mask=a,this._mask&&(this._mask.isMask=!0)}}),Object.defineProperty(d.DisplayObject.prototype,"filters",{get:function(){return this._filters},set:function(a){if(a){for(var b=[],c=0;c<a.length;c++)for(var d=a[c].passes,e=0;e<d.length;e++)b.push(d[e]);this._filterBlock={target:this,filterPasses:b}}this._filters=a}}),d.DisplayObject.prototype.updateTransform=function(){this.rotation!==this.rotationCache&&(this.rotationCache=this.rotation,this._sr=Math.sin(this.rotation),this._cr=Math.cos(this.rotation));var a=this.parent.worldTransform,b=this.worldTransform,c=this.pivot.x,d=this.pivot.y,e=this._cr*this.scale.x,f=-this._sr*this.scale.y,g=this._sr*this.scale.x,h=this._cr*this.scale.y,i=this.position.x-e*c-d*f,j=this.position.y-h*d-c*g,k=a.a,l=a.b,m=a.c,n=a.d;b.a=k*e+l*g,b.b=k*f+l*h,b.tx=k*i+l*j+a.tx,b.c=m*e+n*g,b.d=m*f+n*h,b.ty=m*i+n*j+a.ty,this.worldAlpha=this.alpha*this.parent.worldAlpha},d.DisplayObject.prototype.getBounds=function(a){return a=a,d.EmptyRectangle},d.DisplayObject.prototype.getLocalBounds=function(){return this.getBounds(d.identityMatrix)},d.DisplayObject.prototype.setStageReference=function(a){this.stage=a,this._interactive&&(this.stage.dirty=!0)},d.DisplayObject.prototype._renderWebGL=function(a){a=a},d.DisplayObject.prototype._renderCanvas=function(a){a=a},Object.defineProperty(d.DisplayObject.prototype,"x",{get:function(){return this.position.x},set:function(a){this.position.x=a}}),Object.defineProperty(d.DisplayObject.prototype,"y",{get:function(){return this.position.y},set:function(a){this.position.y=a}}),d.DisplayObjectContainer=function(){d.DisplayObject.call(this),this.children=[]},d.DisplayObjectContainer.prototype=Object.create(d.DisplayObject.prototype),d.DisplayObjectContainer.prototype.constructor=d.DisplayObjectContainer,d.DisplayObjectContainer.prototype.addChild=function(a){this.addChildAt(a,this.children.length)},d.DisplayObjectContainer.prototype.addChildAt=function(a,b){if(!(b>=0&&b<=this.children.length))throw new Error(a+" The index "+b+" supplied is out of bounds "+this.children.length);a.parent&&a.parent.removeChild(a),a.parent=this,this.children.splice(b,0,a),this.stage&&a.setStageReference(this.stage)},d.DisplayObjectContainer.prototype.swapChildren=function(a,b){if(a!==b){var c=this.children.indexOf(a),d=this.children.indexOf(b);if(0>c||0>d)throw new Error("swapChildren: Both the supplied DisplayObjects must be a child of the caller.");this.children[c]=b,this.children[d]=a}},d.DisplayObjectContainer.prototype.getChildAt=function(a){if(a>=0&&a<this.children.length)return this.children[a];throw new Error("The supplied DisplayObjects must be a child of the caller "+this)},d.DisplayObjectContainer.prototype.removeChild=function(a){var b=this.children.indexOf(a);if(-1===b)throw new Error(a+" The supplied DisplayObject must be a child of the caller "+this);this.stage&&a.removeStageReference(),a.parent=void 0,this.children.splice(b,1)},d.DisplayObjectContainer.prototype.updateTransform=function(){if(this.visible){d.DisplayObject.prototype.updateTransform.call(this);for(var a=0,b=this.children.length;b>a;a++)this.children[a].updateTransform()}},d.DisplayObjectContainer.prototype.getBounds=function(a){if(0===this.children.length)return d.EmptyRectangle;if(a){var b=this.worldTransform;this.worldTransform=a,this.updateTransform(),this.worldTransform=b}for(var c,e,f,g=1/0,h=1/0,i=-1/0,j=-1/0,k=!1,l=0,m=this.children.length;m>l;l++){var n=this.children[l];n.visible&&(k=!0,c=this.children[l].getBounds(a),g=g<c.x?g:c.x,h=h<c.y?h:c.y,e=c.width+c.x,f=c.height+c.y,i=i>e?i:e,j=j>f?j:f)}if(!k)return d.EmptyRectangle;var o=this._bounds;return o.x=g,o.y=h,o.width=i-g,o.height=j-h,o},d.DisplayObjectContainer.prototype.getLocalBounds=function(){var a=this.worldTransform;this.worldTransform=d.identityMatrix;for(var b=0,c=this.children.length;c>b;b++)this.children[b].updateTransform();var e=this.getBounds();return this.worldTransform=a,e},d.DisplayObjectContainer.prototype.setStageReference=function(a){this.stage=a,this._interactive&&(this.stage.dirty=!0);for(var b=0,c=this.children.length;c>b;b++){var d=this.children[b];d.setStageReference(a)}},d.DisplayObjectContainer.prototype.removeStageReference=function(){for(var a=0,b=this.children.length;b>a;a++){var c=this.children[a];c.removeStageReference()}this._interactive&&(this.stage.dirty=!0),this.stage=null},d.DisplayObjectContainer.prototype._renderWebGL=function(a){if(this.visible&&!(this.alpha<=0)){var b,c;if(this._mask||this._filters){for(this._mask&&(a.spriteBatch.stop(),a.maskManager.pushMask(this.mask,a),a.spriteBatch.start()),this._filters&&(a.spriteBatch.flush(),a.filterManager.pushFilter(this._filterBlock)),b=0,c=this.children.length;c>b;b++)this.children[b]._renderWebGL(a);a.spriteBatch.stop(),this._filters&&a.filterManager.popFilter(),this._mask&&a.maskManager.popMask(a),a.spriteBatch.start()}else for(b=0,c=this.children.length;c>b;b++)this.children[b]._renderWebGL(a)}},d.DisplayObjectContainer.prototype._renderCanvas=function(a){if(this.visible!==!1&&0!==this.alpha){this._mask&&a.maskManager.pushMask(this._mask,a.context);for(var b=0,c=this.children.length;c>b;b++){var d=this.children[b];d._renderCanvas(a)}this._mask&&a.maskManager.popMask(a.context)}},d.Sprite=function(a){d.DisplayObjectContainer.call(this),this.anchor=new d.Point,this.texture=a,this._width=0,this._height=0,this.tint=16777215,this.blendMode=d.blendModes.NORMAL,a.baseTexture.hasLoaded?this.onTextureUpdate():(this.onTextureUpdateBind=this.onTextureUpdate.bind(this),this.texture.addEventListener("update",this.onTextureUpdateBind)),this.renderable=!0},d.Sprite.prototype=Object.create(d.DisplayObjectContainer.prototype),d.Sprite.prototype.constructor=d.Sprite,Object.defineProperty(d.Sprite.prototype,"width",{get:function(){return this.scale.x*this.texture.frame.width},set:function(a){this.scale.x=a/this.texture.frame.width,this._width=a}}),Object.defineProperty(d.Sprite.prototype,"height",{get:function(){return this.scale.y*this.texture.frame.height},set:function(a){this.scale.y=a/this.texture.frame.height,this._height=a}}),d.Sprite.prototype.setTexture=function(a){this.texture.baseTexture!==a.baseTexture?(this.textureChange=!0,this.texture=a):this.texture=a,this.cachedTint=16777215,this.updateFrame=!0},d.Sprite.prototype.onTextureUpdate=function(){this._width&&(this.scale.x=this._width/this.texture.frame.width),this._height&&(this.scale.y=this._height/this.texture.frame.height),this.updateFrame=!0},d.Sprite.prototype.getBounds=function(a){var b=this.texture.frame.width,c=this.texture.frame.height,d=b*(1-this.anchor.x),e=b*-this.anchor.x,f=c*(1-this.anchor.y),g=c*-this.anchor.y,h=a||this.worldTransform,i=h.a,j=h.c,k=h.b,l=h.d,m=h.tx,n=h.ty,o=i*e+k*g+m,p=l*g+j*e+n,q=i*d+k*g+m,r=l*g+j*d+n,s=i*d+k*f+m,t=l*f+j*d+n,u=i*e+k*f+m,v=l*f+j*e+n,w=-1/0,x=-1/0,y=1/0,z=1/0;y=y>o?o:y,y=y>q?q:y,y=y>s?s:y,y=y>u?u:y,z=z>p?p:z,z=z>r?r:z,z=z>t?t:z,z=z>v?v:z,w=o>w?o:w,w=q>w?q:w,w=s>w?s:w,w=u>w?u:w,x=p>x?p:x,x=r>x?r:x,x=t>x?t:x,x=v>x?v:x;var A=this._bounds;return A.x=y,A.width=w-y,A.y=z,A.height=x-z,this._currentBounds=A,A},d.Sprite.prototype._renderWebGL=function(a){if(this.visible&&!(this.alpha<=0)){var b,c;if(this._mask||this._filters){var d=a.spriteBatch;for(this._mask&&(d.stop(),a.maskManager.pushMask(this.mask,a),d.start()),this._filters&&(d.flush(),a.filterManager.pushFilter(this._filterBlock)),d.render(this),b=0,c=this.children.length;c>b;b++)this.children[b]._renderWebGL(a);d.stop(),this._filters&&a.filterManager.popFilter(),this._mask&&a.maskManager.popMask(a),d.start()}else for(a.spriteBatch.render(this),b=0,c=this.children.length;c>b;b++)this.children[b]._renderWebGL(a)}},d.Sprite.prototype._renderCanvas=function(a){if(this.visible!==!1&&0!==this.alpha){var b=this.texture.frame,c=a.context,e=this.texture;if(this.blendMode!==a.currentBlendMode&&(a.currentBlendMode=this.blendMode,c.globalCompositeOperation=d.blendModesCanvas[a.currentBlendMode]),this._mask&&a.maskManager.pushMask(this._mask,a.context),b&&b.width&&b.height&&e.baseTexture.source){c.globalAlpha=this.worldAlpha;var f=this.worldTransform;if(a.roundPixels?c.setTransform(f.a,f.c,f.b,f.d,f.tx||0,f.ty||0):c.setTransform(f.a,f.c,f.b,f.d,f.tx,f.ty),a.smoothProperty&&a.scaleMode!==this.texture.baseTexture.scaleMode&&(a.scaleMode=this.texture.baseTexture.scaleMode,c[a.smoothProperty]=a.scaleMode===d.scaleModes.LINEAR),16777215!==this.tint){if(this.cachedTint!==this.tint){if(!e.baseTexture.hasLoaded)return;this.cachedTint=this.tint,this.tintedTexture=d.CanvasTinter.getTintedTexture(this,this.tint)}c.drawImage(this.tintedTexture,0,0,b.width,b.height,this.anchor.x*-b.width,this.anchor.y*-b.height,b.width,b.height)}else if(e.trim){var g=e.trim;c.drawImage(this.texture.baseTexture.source,b.x,b.y,b.width,b.height,g.x-this.anchor.x*g.width,g.y-this.anchor.y*g.height,b.width,b.height)}else c.drawImage(this.texture.baseTexture.source,b.x,b.y,b.width,b.height,this.anchor.x*-b.width,this.anchor.y*-b.height,b.width,b.height)}for(var h=0,i=this.children.length;i>h;h++){var j=this.children[h];j._renderCanvas(a)}this._mask&&a.maskManager.popMask(a.context)}},d.Sprite.fromFrame=function(a){var b=d.TextureCache[a];if(!b)throw new Error('The frameId "'+a+'" does not exist in the texture cache'+this);return new d.Sprite(b)},d.Sprite.fromImage=function(a,b,c){var e=d.Texture.fromImage(a,b,c);return new d.Sprite(e)},d.SpriteBatch=function(a){d.DisplayObjectContainer.call(this),this.textureThing=a,this.ready=!1},d.SpriteBatch.prototype=Object.create(d.DisplayObjectContainer.prototype),d.SpriteBatch.constructor=d.SpriteBatch,d.SpriteBatch.prototype.initWebGL=function(a){this.fastSpriteBatch=new d.WebGLFastSpriteBatch(a),this.ready=!0},d.SpriteBatch.prototype.updateTransform=function(){d.DisplayObject.prototype.updateTransform.call(this)},d.SpriteBatch.prototype._renderWebGL=function(a){!this.visible||this.alpha<=0||!this.children.length||(this.ready||this.initWebGL(a.gl),a.spriteBatch.stop(),a.shaderManager.activateShader(a.shaderManager.fastShader),this.fastSpriteBatch.begin(this,a),this.fastSpriteBatch.render(this),a.shaderManager.activateShader(a.shaderManager.defaultShader),a.spriteBatch.start())},d.SpriteBatch.prototype._renderCanvas=function(a){var b=a.context;b.globalAlpha=this.worldAlpha;var c=this.worldTransform;a.roundPixels?b.setTransform(c.a,c.c,c.b,c.d,Math.floor(c.tx),Math.floor(c.ty)):b.setTransform(c.a,c.c,c.b,c.d,c.tx,c.ty),b.save();for(var e=0;e<this.children.length;e++){var f=this.children[e],g=f.texture,h=g.frame;if(b.globalAlpha=this.worldAlpha*f.alpha,f.rotation%(2*Math.PI)===0)b.drawImage(g.baseTexture.source,h.x,h.y,h.width,h.height,f.anchor.x*-h.width*f.scale.x+f.position.x+.5|0,f.anchor.y*-h.height*f.scale.y+f.position.y+.5|0,h.width*f.scale.x,h.height*f.scale.y);else{d.DisplayObject.prototype.updateTransform.call(f),c=f.localTransform,this.rotation!==this.rotationCache&&(this.rotationCache=this.rotation,this._sr=Math.sin(this.rotation),this._cr=Math.cos(this.rotation));var i=f._cr*f.scale.x,j=-f._sr*f.scale.y,k=f._sr*f.scale.x,l=f._cr*f.scale.y;b.setTransform(i,k,j,l,f.position.x,f.position.y),b.drawImage(g.baseTexture.source,h.x,h.y,h.width,h.height,f.anchor.x*-h.width+.5|0,f.anchor.y*-h.height+.5|0,h.width,h.height)}}b.restore()},d.MovieClip=function(a){d.Sprite.call(this,a[0]),this.textures=a,this.animationSpeed=1,this.loop=!0,this.onComplete=null,this.currentFrame=0,this.playing=!1},d.MovieClip.prototype=Object.create(d.Sprite.prototype),d.MovieClip.prototype.constructor=d.MovieClip,Object.defineProperty(d.MovieClip.prototype,"totalFrames",{get:function(){return this.textures.length}}),d.MovieClip.prototype.stop=function(){this.playing=!1},d.MovieClip.prototype.play=function(){this.playing=!0},d.MovieClip.prototype.gotoAndStop=function(a){this.playing=!1,this.currentFrame=a;var b=this.currentFrame+.5|0;this.setTexture(this.textures[b%this.textures.length])},d.MovieClip.prototype.gotoAndPlay=function(a){this.currentFrame=a,this.playing=!0},d.MovieClip.prototype.updateTransform=function(){if(d.Sprite.prototype.updateTransform.call(this),this.playing){this.currentFrame+=this.animationSpeed;var a=this.currentFrame+.5|0;this.loop||a<this.textures.length?this.setTexture(this.textures[a%this.textures.length]):a>=this.textures.length&&(this.gotoAndStop(this.textures.length-1),this.onComplete&&this.onComplete())}},d.FilterBlock=function(){this.visible=!0,this.renderable=!0},d.Text=function(a,b){this.canvas=document.createElement("canvas"),this.context=this.canvas.getContext("2d"),d.Sprite.call(this,d.Texture.fromCanvas(this.canvas)),this.setText(a),this.setStyle(b),this.updateText(),this.dirty=!1},d.Text.prototype=Object.create(d.Sprite.prototype),d.Text.prototype.constructor=d.Text,d.Text.prototype.setStyle=function(a){a=a||{},a.font=a.font||"bold 20pt Arial",a.fill=a.fill||"black",a.align=a.align||"left",a.stroke=a.stroke||"black",a.strokeThickness=a.strokeThickness||0,a.wordWrap=a.wordWrap||!1,a.wordWrapWidth=a.wordWrapWidth||100,this.style=a,this.dirty=!0},d.Text.prototype.setText=function(a){this.text=a.toString()||" ",this.dirty=!0},d.Text.prototype.updateText=function(){this.context.font=this.style.font;var a=this.text;this.style.wordWrap&&(a=this.wordWrap(this.text));for(var b=a.split(/(?:\r\n|\r|\n)/),c=[],e=0,f=0;f<b.length;f++){var g=this.context.measureText(b[f]).width;c[f]=g,e=Math.max(e,g)}this.canvas.width=e+this.style.strokeThickness;var h=this.determineFontHeight("font: "+this.style.font+";")+this.style.strokeThickness;for(this.canvas.height=h*b.length,navigator.isCocoonJS&&this.context.clearRect(0,0,this.canvas.width,this.canvas.height),this.context.fillStyle=this.style.fill,this.context.font=this.style.font,this.context.strokeStyle=this.style.stroke,this.context.lineWidth=this.style.strokeThickness,this.context.textBaseline="top",f=0;f<b.length;f++){var i=new d.Point(this.style.strokeThickness/2,this.style.strokeThickness/2+f*h);"right"===this.style.align?i.x+=e-c[f]:"center"===this.style.align&&(i.x+=(e-c[f])/2),this.style.stroke&&this.style.strokeThickness&&this.context.strokeText(b[f],i.x,i.y),this.style.fill&&this.context.fillText(b[f],i.x,i.y)}this.updateTexture()},d.Text.prototype.updateTexture=function(){this.texture.baseTexture.width=this.canvas.width,this.texture.baseTexture.height=this.canvas.height,this.texture.frame.width=this.canvas.width,this.texture.frame.height=this.canvas.height,this._width=this.canvas.width,this._height=this.canvas.height,this.requiresUpdate=!0},d.Text.prototype._renderWebGL=function(a){this.requiresUpdate&&(this.requiresUpdate=!1,d.updateWebGLTexture(this.texture.baseTexture,a.gl)),d.Sprite.prototype._renderWebGL.call(this,a)},d.Text.prototype.updateTransform=function(){this.dirty&&(this.updateText(),this.dirty=!1),d.Sprite.prototype.updateTransform.call(this)},d.Text.prototype.determineFontHeight=function(a){var b=d.Text.heightCache[a];if(!b){var c=document.getElementsByTagName("body")[0],e=document.createElement("div"),f=document.createTextNode("M");e.appendChild(f),e.setAttribute("style",a+";position:absolute;top:0;left:0"),c.appendChild(e),b=e.offsetHeight,d.Text.heightCache[a]=b,c.removeChild(e)}return b},d.Text.prototype.wordWrap=function(a){for(var b="",c=a.split("\n"),d=0;d<c.length;d++){for(var e=this.style.wordWrapWidth,f=c[d].split(" "),g=0;g<f.length;g++){var h=this.context.measureText(f[g]).width,i=h+this.context.measureText(" ").width;i>e?(g>0&&(b+="\n"),b+=f[g]+" ",e=this.style.wordWrapWidth-h):(e-=i,b+=f[g]+" ")}d<c.length-1&&(b+="\n")}return b},d.Text.prototype.destroy=function(a){a&&this.texture.destroy()},d.Text.heightCache={},d.BitmapText=function(a,b){d.DisplayObjectContainer.call(this),this._pool=[],this.setText(a),this.setStyle(b),this.updateText(),this.dirty=!1},d.BitmapText.prototype=Object.create(d.DisplayObjectContainer.prototype),d.BitmapText.prototype.constructor=d.BitmapText,d.BitmapText.prototype.setText=function(a){this.text=a||" ",this.dirty=!0},d.BitmapText.prototype.setStyle=function(a){a=a||{},a.align=a.align||"left",this.style=a;var b=a.font.split(" ");this.fontName=b[b.length-1],this.fontSize=b.length>=2?parseInt(b[b.length-2],10):d.BitmapText.fonts[this.fontName].size,this.dirty=!0,this.tint=a.tint},d.BitmapText.prototype.updateText=function(){for(var a=d.BitmapText.fonts[this.fontName],b=new d.Point,c=null,e=[],f=0,g=[],h=0,i=this.fontSize/a.size,j=0;j<this.text.length;j++){var k=this.text.charCodeAt(j);if(/(?:\r\n|\r|\n)/.test(this.text.charAt(j)))g.push(b.x),f=Math.max(f,b.x),h++,b.x=0,b.y+=a.lineHeight,c=null;else{var l=a.chars[k];l&&(c&&l[c]&&(b.x+=l.kerning[c]),e.push({texture:l.texture,line:h,charCode:k,position:new d.Point(b.x+l.xOffset,b.y+l.yOffset)}),b.x+=l.xAdvance,c=k)}}g.push(b.x),f=Math.max(f,b.x);var m=[];for(j=0;h>=j;j++){var n=0;"right"===this.style.align?n=f-g[j]:"center"===this.style.align&&(n=(f-g[j])/2),m.push(n)}var o=this.children.length,p=e.length,q=this.tint||16777215;for(j=0;p>j;j++){var r=o>j?this.children[j]:this._pool.pop();r?r.setTexture(e[j].texture):r=new d.Sprite(e[j].texture),r.position.x=(e[j].position.x+m[e[j].line])*i,r.position.y=e[j].position.y*i,r.scale.x=r.scale.y=i,r.tint=q,r.parent||this.addChild(r)}for(;this.children.length>p;){var s=this.getChildAt(this.children.length-1);this._pool.push(s),this.removeChild(s)}this.textWidth=f*i,this.textHeight=(b.y+a.lineHeight)*i},d.BitmapText.prototype.updateTransform=function(){this.dirty&&(this.updateText(),this.dirty=!1),d.DisplayObjectContainer.prototype.updateTransform.call(this)},d.BitmapText.fonts={},d.InteractionData=function(){this.global=new d.Point,this.local=new d.Point,this.target=null,this.originalEvent=null},d.InteractionData.prototype.getLocalPosition=function(a){var b=a.worldTransform,c=this.global,e=b.a,f=b.b,g=b.tx,h=b.c,i=b.d,j=b.ty,k=1/(e*i+f*-h);return new d.Point(i*k*c.x+-f*k*c.y+(j*f-g*i)*k,e*k*c.y+-h*k*c.x+(-j*e+g*h)*k)},d.InteractionData.prototype.constructor=d.InteractionData,d.InteractionManager=function(a){this.stage=a,this.mouse=new d.InteractionData,this.touchs={},this.tempPoint=new d.Point,this.mouseoverEnabled=!0,this.pool=[],this.interactiveItems=[],this.interactionDOMElement=null,this.onMouseMove=this.onMouseMove.bind(this),this.onMouseDown=this.onMouseDown.bind(this),this.onMouseOut=this.onMouseOut.bind(this),this.onMouseUp=this.onMouseUp.bind(this),this.onTouchStart=this.onTouchStart.bind(this),this.onTouchEnd=this.onTouchEnd.bind(this),this.onTouchMove=this.onTouchMove.bind(this),this.last=0,this.currentCursorStyle="inherit",this.mouseOut=!1},d.InteractionManager.prototype.constructor=d.InteractionManager,d.InteractionManager.prototype.collectInteractiveSprite=function(a,b){for(var c=a.children,d=c.length,e=d-1;e>=0;e--){var f=c[e];f.interactive?(b.interactiveChildren=!0,this.interactiveItems.push(f),f.children.length>0&&this.collectInteractiveSprite(f,f)):(f.__iParent=null,f.children.length>0&&this.collectInteractiveSprite(f,b))}},d.InteractionManager.prototype.setTarget=function(a){this.target=a,null===this.interactionDOMElement&&this.setTargetDomElement(a.view)},d.InteractionManager.prototype.setTargetDomElement=function(a){this.removeEvents(),window.navigator.msPointerEnabled&&(a.style["-ms-content-zooming"]="none",a.style["-ms-touch-action"]="none"),this.interactionDOMElement=a,a.addEventListener("mousemove",this.onMouseMove,!0),a.addEventListener("mousedown",this.onMouseDown,!0),a.addEventListener("mouseout",this.onMouseOut,!0),a.addEventListener("touchstart",this.onTouchStart,!0),a.addEventListener("touchend",this.onTouchEnd,!0),a.addEventListener("touchmove",this.onTouchMove,!0),document.body.addEventListener("mouseup",this.onMouseUp,!0)},d.InteractionManager.prototype.removeEvents=function(){this.interactionDOMElement&&(this.interactionDOMElement.style["-ms-content-zooming"]="",this.interactionDOMElement.style["-ms-touch-action"]="",this.interactionDOMElement.removeEventListener("mousemove",this.onMouseMove,!0),this.interactionDOMElement.removeEventListener("mousedown",this.onMouseDown,!0),this.interactionDOMElement.removeEventListener("mouseout",this.onMouseOut,!0),this.interactionDOMElement.removeEventListener("touchstart",this.onTouchStart,!0),this.interactionDOMElement.removeEventListener("touchend",this.onTouchEnd,!0),this.interactionDOMElement.removeEventListener("touchmove",this.onTouchMove,!0),this.interactionDOMElement=null,document.body.removeEventListener("mouseup",this.onMouseUp,!0))},d.InteractionManager.prototype.update=function(){if(this.target){var a=Date.now(),b=a-this.last;if(b=b*d.INTERACTION_FREQUENCY/1e3,!(1>b)){this.last=a;var c=0;if(this.dirty){this.dirty=!1;var e=this.interactiveItems.length;for(c=0;e>c;c++)this.interactiveItems[c].interactiveChildren=!1;this.interactiveItems=[],this.stage.interactive&&this.interactiveItems.push(this.stage),this.collectInteractiveSprite(this.stage,this.stage)}var f=this.interactiveItems.length,g="inherit",h=!1;for(c=0;f>c;c++){var i=this.interactiveItems[c];i.__hit=this.hitTest(i,this.mouse),this.mouse.target=i,i.__hit&&!h?(i.buttonMode&&(g=i.defaultCursor),i.interactiveChildren||(h=!0),i.__isOver||(i.mouseover&&i.mouseover(this.mouse),i.__isOver=!0)):i.__isOver&&(i.mouseout&&i.mouseout(this.mouse),i.__isOver=!1)}this.currentCursorStyle!==g&&(this.currentCursorStyle=g,this.interactionDOMElement.style.cursor=g)}}},d.InteractionManager.prototype.onMouseMove=function(a){this.mouse.originalEvent=a||window.event;var b=this.interactionDOMElement.getBoundingClientRect();this.mouse.global.x=(a.clientX-b.left)*(this.target.width/b.width),this.mouse.global.y=(a.clientY-b.top)*(this.target.height/b.height);for(var c=this.interactiveItems.length,d=0;c>d;d++){var e=this.interactiveItems[d];e.mousemove&&e.mousemove(this.mouse)}},d.InteractionManager.prototype.onMouseDown=function(a){this.mouse.originalEvent=a||window.event,d.AUTO_PREVENT_DEFAULT&&this.mouse.originalEvent.preventDefault();for(var b=this.interactiveItems.length,c=0;b>c;c++){var e=this.interactiveItems[c];if((e.mousedown||e.click)&&(e.__mouseIsDown=!0,e.__hit=this.hitTest(e,this.mouse),e.__hit&&(e.mousedown&&e.mousedown(this.mouse),e.__isDown=!0,!e.interactiveChildren)))break}},d.InteractionManager.prototype.onMouseOut=function(){var a=this.interactiveItems.length;this.interactionDOMElement.style.cursor="inherit";for(var b=0;a>b;b++){var c=this.interactiveItems[b];c.__isOver&&(this.mouse.target=c,c.mouseout&&c.mouseout(this.mouse),c.__isOver=!1)}this.mouseOut=!0,this.mouse.global.x=-1e4,this.mouse.global.y=-1e4},d.InteractionManager.prototype.onMouseUp=function(a){this.mouse.originalEvent=a||window.event;for(var b=this.interactiveItems.length,c=!1,d=0;b>d;d++){var e=this.interactiveItems[d];e.__hit=this.hitTest(e,this.mouse),e.__hit&&!c?(e.mouseup&&e.mouseup(this.mouse),e.__isDown&&e.click&&e.click(this.mouse),e.interactiveChildren||(c=!0)):e.__isDown&&e.mouseupoutside&&e.mouseupoutside(this.mouse),e.__isDown=!1}},d.InteractionManager.prototype.hitTest=function(a,b){var c=b.global;if(!a.worldVisible)return!1;var e=a instanceof d.Sprite,f=a.worldTransform,g=f.a,h=f.b,i=f.tx,j=f.c,k=f.d,l=f.ty,m=1/(g*k+h*-j),n=k*m*c.x+-h*m*c.y+(l*h-i*k)*m,o=g*m*c.y+-j*m*c.x+(-l*g+i*j)*m;if(b.target=a,a.hitArea&&a.hitArea.contains)return a.hitArea.contains(n,o)?(b.target=a,!0):!1;if(e){var p,q=a.texture.frame.width,r=a.texture.frame.height,s=-q*a.anchor.x;if(n>s&&s+q>n&&(p=-r*a.anchor.y,o>p&&p+r>o))return b.target=a,!0}for(var t=a.children.length,u=0;t>u;u++){var v=a.children[u],w=this.hitTest(v,b);if(w)return b.target=a,!0}return!1},d.InteractionManager.prototype.onTouchMove=function(a){var b,c=this.interactionDOMElement.getBoundingClientRect(),d=a.changedTouches,e=0;for(e=0;e<d.length;e++){var f=d[e];b=this.touchs[f.identifier],b.originalEvent=a||window.event,b.global.x=(f.clientX-c.left)*(this.target.width/c.width),b.global.y=(f.clientY-c.top)*(this.target.height/c.height),navigator.isCocoonJS&&(b.global.x=f.clientX,b.global.y=f.clientY)}var g=this.interactiveItems.length;for(e=0;g>e;e++){var h=this.interactiveItems[e];h.touchmove&&h.touchmove(b)}},d.InteractionManager.prototype.onTouchStart=function(a){var b=this.interactionDOMElement.getBoundingClientRect();d.AUTO_PREVENT_DEFAULT&&a.preventDefault();for(var c=a.changedTouches,e=0;e<c.length;e++){var f=c[e],g=this.pool.pop();g||(g=new d.InteractionData),g.originalEvent=a||window.event,this.touchs[f.identifier]=g,g.global.x=(f.clientX-b.left)*(this.target.width/b.width),g.global.y=(f.clientY-b.top)*(this.target.height/b.height),navigator.isCocoonJS&&(g.global.x=f.clientX,g.global.y=f.clientY);for(var h=this.interactiveItems.length,i=0;h>i;i++){var j=this.interactiveItems[i];if((j.touchstart||j.tap)&&(j.__hit=this.hitTest(j,g),j.__hit&&(j.touchstart&&j.touchstart(g),j.__isDown=!0,j.__touchData=g,!j.interactiveChildren)))break}}},d.InteractionManager.prototype.onTouchEnd=function(a){for(var b=this.interactionDOMElement.getBoundingClientRect(),c=a.changedTouches,d=0;d<c.length;d++){var e=c[d],f=this.touchs[e.identifier],g=!1;f.global.x=(e.clientX-b.left)*(this.target.width/b.width),f.global.y=(e.clientY-b.top)*(this.target.height/b.height),navigator.isCocoonJS&&(f.global.x=e.clientX,f.global.y=e.clientY);for(var h=this.interactiveItems.length,i=0;h>i;i++){var j=this.interactiveItems[i],k=j.__touchData;j.__hit=this.hitTest(j,f),k===f&&(f.originalEvent=a||window.event,(j.touchend||j.tap)&&(j.__hit&&!g?(j.touchend&&j.touchend(f),j.__isDown&&j.tap&&j.tap(f),j.interactiveChildren||(g=!0)):j.__isDown&&j.touchendoutside&&j.touchendoutside(f),j.__isDown=!1),j.__touchData=null)}this.pool.push(f),this.touchs[e.identifier]=null}},d.Stage=function(a){d.DisplayObjectContainer.call(this),this.worldTransform=new d.Matrix,this.interactive=!0,this.interactionManager=new d.InteractionManager(this),this.dirty=!0,this.stage=this,this.stage.hitArea=new d.Rectangle(0,0,1e5,1e5),this.setBackgroundColor(a)},d.Stage.prototype=Object.create(d.DisplayObjectContainer.prototype),d.Stage.prototype.constructor=d.Stage,d.Stage.prototype.setInteractionDelegate=function(a){this.interactionManager.setTargetDomElement(a)},d.Stage.prototype.updateTransform=function(){this.worldAlpha=1;for(var a=0,b=this.children.length;b>a;a++)this.children[a].updateTransform();this.dirty&&(this.dirty=!1,this.interactionManager.dirty=!0),this.interactive&&this.interactionManager.update()},d.Stage.prototype.setBackgroundColor=function(a){this.backgroundColor=a||0,this.backgroundColorSplit=d.hex2rgb(this.backgroundColor);var b=this.backgroundColor.toString(16);b="000000".substr(0,6-b.length)+b,this.backgroundColorString="#"+b},d.Stage.prototype.getMousePosition=function(){return this.interactionManager.mouse.global};for(var e=0,f=["ms","moz","webkit","o"],h=0;h<f.length&&!window.requestAnimationFrame;++h)window.requestAnimationFrame=window[f[h]+"RequestAnimationFrame"],window.cancelAnimationFrame=window[f[h]+"CancelAnimationFrame"]||window[f[h]+"CancelRequestAnimationFrame"];window.requestAnimationFrame||(window.requestAnimationFrame=function(a){var b=(new Date).getTime(),c=Math.max(0,16-(b-e)),d=window.setTimeout(function(){a(b+c)},c);return e=b+c,d}),window.cancelAnimationFrame||(window.cancelAnimationFrame=function(a){clearTimeout(a)}),window.requestAnimFrame=window.requestAnimationFrame,d.hex2rgb=function(a){return[(a>>16&255)/255,(a>>8&255)/255,(255&a)/255]},d.rgb2hex=function(a){return(255*a[0]<<16)+(255*a[1]<<8)+255*a[2]},"function"!=typeof Function.prototype.bind&&(Function.prototype.bind=function(){var a=Array.prototype.slice;return function(b){function c(){var f=e.concat(a.call(arguments));
d.apply(this instanceof c?this:b,f)}var d=this,e=a.call(arguments,1);if("function"!=typeof d)throw new TypeError;return c.prototype=function f(a){return a&&(f.prototype=a),this instanceof f?void 0:new f}(d.prototype),c}}()),d.AjaxRequest=function(){var a=["Msxml2.XMLHTTP.6.0","Msxml2.XMLHTTP.3.0","Microsoft.XMLHTTP"];if(!window.ActiveXObject)return window.XMLHttpRequest?new window.XMLHttpRequest:!1;for(var b=0;b<a.length;b++)try{return new window.ActiveXObject(a[b])}catch(c){}},d.canUseNewCanvasBlendModes=function(){var a=document.createElement("canvas");a.width=1,a.height=1;var b=a.getContext("2d");return b.fillStyle="#000",b.fillRect(0,0,1,1),b.globalCompositeOperation="multiply",b.fillStyle="#fff",b.fillRect(0,0,1,1),0===b.getImageData(0,0,1,1).data[0]},d.getNextPowerOfTwo=function(a){if(a>0&&0===(a&a-1))return a;for(var b=1;a>b;)b<<=1;return b},d.EventTarget=function(){var a={};this.addEventListener=this.on=function(b,c){void 0===a[b]&&(a[b]=[]),-1===a[b].indexOf(c)&&a[b].push(c)},this.dispatchEvent=this.emit=function(b){if(a[b.type]&&a[b.type].length)for(var c=0,d=a[b.type].length;d>c;c++)a[b.type][c](b)},this.removeEventListener=this.off=function(b,c){var d=a[b].indexOf(c);-1!==d&&a[b].splice(d,1)},this.removeAllEventListeners=function(b){var c=a[b];c&&(c.length=0)}},d.autoDetectRenderer=function(a,b,c,e,f){a||(a=800),b||(b=600);var g=function(){try{var a=document.createElement("canvas");return!!window.WebGLRenderingContext&&(a.getContext("webgl")||a.getContext("experimental-webgl"))}catch(b){return!1}}();return g?new d.WebGLRenderer(a,b,c,e,f):new d.CanvasRenderer(a,b,c,e)},d.PolyK={},d.PolyK.Triangulate=function(a){var b=!0,c=a.length>>1;if(3>c)return[];for(var e=[],f=[],g=0;c>g;g++)f.push(g);g=0;for(var h=c;h>3;){var i=f[(g+0)%h],j=f[(g+1)%h],k=f[(g+2)%h],l=a[2*i],m=a[2*i+1],n=a[2*j],o=a[2*j+1],p=a[2*k],q=a[2*k+1],r=!1;if(d.PolyK._convex(l,m,n,o,p,q,b)){r=!0;for(var s=0;h>s;s++){var t=f[s];if(t!==i&&t!==j&&t!==k&&d.PolyK._PointInTriangle(a[2*t],a[2*t+1],l,m,n,o,p,q)){r=!1;break}}}if(r)e.push(i,j,k),f.splice((g+1)%h,1),h--,g=0;else if(g++>3*h){if(!b)return window.console.log("PIXI Warning: shape too complex to fill"),[];for(e=[],f=[],g=0;c>g;g++)f.push(g);g=0,h=c,b=!1}}return e.push(f[0],f[1],f[2]),e},d.PolyK._PointInTriangle=function(a,b,c,d,e,f,g,h){var i=g-c,j=h-d,k=e-c,l=f-d,m=a-c,n=b-d,o=i*i+j*j,p=i*k+j*l,q=i*m+j*n,r=k*k+l*l,s=k*m+l*n,t=1/(o*r-p*p),u=(r*q-p*s)*t,v=(o*s-p*q)*t;return u>=0&&v>=0&&1>u+v},d.PolyK._convex=function(a,b,c,d,e,f,g){return(b-d)*(e-c)+(c-a)*(f-d)>=0===g},d.initDefaultShaders=function(){},d.CompileVertexShader=function(a,b){return d._CompileShader(a,b,a.VERTEX_SHADER)},d.CompileFragmentShader=function(a,b){return d._CompileShader(a,b,a.FRAGMENT_SHADER)},d._CompileShader=function(a,b,c){var d=b.join("\n"),e=a.createShader(c);return a.shaderSource(e,d),a.compileShader(e),a.getShaderParameter(e,a.COMPILE_STATUS)?e:(window.console.log(a.getShaderInfoLog(e)),null)},d.compileProgram=function(a,b,c){var e=d.CompileFragmentShader(a,c),f=d.CompileVertexShader(a,b),g=a.createProgram();return a.attachShader(g,f),a.attachShader(g,e),a.linkProgram(g),a.getProgramParameter(g,a.LINK_STATUS)||window.console.log("Could not initialise shaders"),g},d.PixiShader=function(a){this.gl=a,this.program=null,this.fragmentSrc=["precision lowp float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform sampler2D uSampler;","void main(void) {","   gl_FragColor = texture2D(uSampler, vTextureCoord) * vColor ;","}"],this.textureCount=0,this.attributes=[],this.init()},d.PixiShader.prototype.init=function(){var a=this.gl,b=d.compileProgram(a,this.vertexSrc||d.PixiShader.defaultVertexSrc,this.fragmentSrc);a.useProgram(b),this.uSampler=a.getUniformLocation(b,"uSampler"),this.projectionVector=a.getUniformLocation(b,"projectionVector"),this.offsetVector=a.getUniformLocation(b,"offsetVector"),this.dimensions=a.getUniformLocation(b,"dimensions"),this.aVertexPosition=a.getAttribLocation(b,"aVertexPosition"),this.aTextureCoord=a.getAttribLocation(b,"aTextureCoord"),this.colorAttribute=a.getAttribLocation(b,"aColor"),-1===this.colorAttribute&&(this.colorAttribute=2),this.attributes=[this.aVertexPosition,this.aTextureCoord,this.colorAttribute];for(var c in this.uniforms)this.uniforms[c].uniformLocation=a.getUniformLocation(b,c);this.initUniforms(),this.program=b},d.PixiShader.prototype.initUniforms=function(){this.textureCount=1;var a,b=this.gl;for(var c in this.uniforms){a=this.uniforms[c];var d=a.type;"sampler2D"===d?(a._init=!1,null!==a.value&&this.initSampler2D(a)):"mat2"===d||"mat3"===d||"mat4"===d?(a.glMatrix=!0,a.glValueLength=1,"mat2"===d?a.glFunc=b.uniformMatrix2fv:"mat3"===d?a.glFunc=b.uniformMatrix3fv:"mat4"===d&&(a.glFunc=b.uniformMatrix4fv)):(a.glFunc=b["uniform"+d],a.glValueLength="2f"===d||"2i"===d?2:"3f"===d||"3i"===d?3:"4f"===d||"4i"===d?4:1)}},d.PixiShader.prototype.initSampler2D=function(a){if(a.value&&a.value.baseTexture&&a.value.baseTexture.hasLoaded){var b=this.gl;if(b.activeTexture(b["TEXTURE"+this.textureCount]),b.bindTexture(b.TEXTURE_2D,a.value.baseTexture._glTexture),a.textureData){var c=a.textureData,d=c.magFilter?c.magFilter:b.LINEAR,e=c.minFilter?c.minFilter:b.LINEAR,f=c.wrapS?c.wrapS:b.CLAMP_TO_EDGE,g=c.wrapT?c.wrapT:b.CLAMP_TO_EDGE,h=c.luminance?b.LUMINANCE:b.RGBA;if(c.repeat&&(f=b.REPEAT,g=b.REPEAT),b.pixelStorei(b.UNPACK_FLIP_Y_WEBGL,!!c.flipY),c.width){var i=c.width?c.width:512,j=c.height?c.height:2,k=c.border?c.border:0;b.texImage2D(b.TEXTURE_2D,0,h,i,j,k,h,b.UNSIGNED_BYTE,null)}else b.texImage2D(b.TEXTURE_2D,0,h,b.RGBA,b.UNSIGNED_BYTE,a.value.baseTexture.source);b.texParameteri(b.TEXTURE_2D,b.TEXTURE_MAG_FILTER,d),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_MIN_FILTER,e),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_S,f),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_T,g)}b.uniform1i(a.uniformLocation,this.textureCount),a._init=!0,this.textureCount++}},d.PixiShader.prototype.syncUniforms=function(){this.textureCount=1;var a,b=this.gl;for(var c in this.uniforms)a=this.uniforms[c],1===a.glValueLength?a.glMatrix===!0?a.glFunc.call(b,a.uniformLocation,a.transpose,a.value):a.glFunc.call(b,a.uniformLocation,a.value):2===a.glValueLength?a.glFunc.call(b,a.uniformLocation,a.value.x,a.value.y):3===a.glValueLength?a.glFunc.call(b,a.uniformLocation,a.value.x,a.value.y,a.value.z):4===a.glValueLength?a.glFunc.call(b,a.uniformLocation,a.value.x,a.value.y,a.value.z,a.value.w):"sampler2D"===a.type&&(a._init?(b.activeTexture(b["TEXTURE"+this.textureCount]),b.bindTexture(b.TEXTURE_2D,a.value.baseTexture._glTextures[b.id]||d.createWebGLTexture(a.value.baseTexture,b)),b.uniform1i(a.uniformLocation,this.textureCount),this.textureCount++):this.initSampler2D(a))},d.PixiShader.prototype.destroy=function(){this.gl.deleteProgram(this.program),this.uniforms=null,this.gl=null,this.attributes=null},d.PixiShader.defaultVertexSrc=["attribute vec2 aVertexPosition;","attribute vec2 aTextureCoord;","attribute vec2 aColor;","uniform vec2 projectionVector;","uniform vec2 offsetVector;","varying vec2 vTextureCoord;","varying vec4 vColor;","const vec2 center = vec2(-1.0, 1.0);","void main(void) {","   gl_Position = vec4( ((aVertexPosition + offsetVector) / projectionVector) + center , 0.0, 1.0);","   vTextureCoord = aTextureCoord;","   vec3 color = mod(vec3(aColor.y/65536.0, aColor.y/256.0, aColor.y), 256.0) / 256.0;","   vColor = vec4(color * aColor.x, aColor.x);","}"],d.PixiFastShader=function(a){this.gl=a,this.program=null,this.fragmentSrc=["precision lowp float;","varying vec2 vTextureCoord;","varying float vColor;","uniform sampler2D uSampler;","void main(void) {","   gl_FragColor = texture2D(uSampler, vTextureCoord) * vColor ;","}"],this.vertexSrc=["attribute vec2 aVertexPosition;","attribute vec2 aPositionCoord;","attribute vec2 aScale;","attribute float aRotation;","attribute vec2 aTextureCoord;","attribute float aColor;","uniform vec2 projectionVector;","uniform vec2 offsetVector;","uniform mat3 uMatrix;","varying vec2 vTextureCoord;","varying float vColor;","const vec2 center = vec2(-1.0, 1.0);","void main(void) {","   vec2 v;","   vec2 sv = aVertexPosition * aScale;","   v.x = (sv.x) * cos(aRotation) - (sv.y) * sin(aRotation);","   v.y = (sv.x) * sin(aRotation) + (sv.y) * cos(aRotation);","   v = ( uMatrix * vec3(v + aPositionCoord , 1.0) ).xy ;","   gl_Position = vec4( ( v / projectionVector) + center , 0.0, 1.0);","   vTextureCoord = aTextureCoord;","   vColor = aColor;","}"],this.textureCount=0,this.init()},d.PixiFastShader.prototype.init=function(){var a=this.gl,b=d.compileProgram(a,this.vertexSrc,this.fragmentSrc);a.useProgram(b),this.uSampler=a.getUniformLocation(b,"uSampler"),this.projectionVector=a.getUniformLocation(b,"projectionVector"),this.offsetVector=a.getUniformLocation(b,"offsetVector"),this.dimensions=a.getUniformLocation(b,"dimensions"),this.uMatrix=a.getUniformLocation(b,"uMatrix"),this.aVertexPosition=a.getAttribLocation(b,"aVertexPosition"),this.aPositionCoord=a.getAttribLocation(b,"aPositionCoord"),this.aScale=a.getAttribLocation(b,"aScale"),this.aRotation=a.getAttribLocation(b,"aRotation"),this.aTextureCoord=a.getAttribLocation(b,"aTextureCoord"),this.colorAttribute=a.getAttribLocation(b,"aColor"),-1===this.colorAttribute&&(this.colorAttribute=2),this.attributes=[this.aVertexPosition,this.aPositionCoord,this.aScale,this.aRotation,this.aTextureCoord,this.colorAttribute],this.program=b},d.PixiFastShader.prototype.destroy=function(){this.gl.deleteProgram(this.program),this.uniforms=null,this.gl=null,this.attributes=null},d.StripShader=function(){this.program=null,this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying float vColor;","uniform float alpha;","uniform sampler2D uSampler;","void main(void) {","   gl_FragColor = texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y));","   gl_FragColor = gl_FragColor * alpha;","}"],this.vertexSrc=["attribute vec2 aVertexPosition;","attribute vec2 aTextureCoord;","attribute float aColor;","uniform mat3 translationMatrix;","uniform vec2 projectionVector;","varying vec2 vTextureCoord;","uniform vec2 offsetVector;","varying float vColor;","void main(void) {","   vec3 v = translationMatrix * vec3(aVertexPosition, 1.0);","   v -= offsetVector.xyx;","   gl_Position = vec4( v.x / projectionVector.x -1.0, v.y / projectionVector.y + 1.0 , 0.0, 1.0);","   vTextureCoord = aTextureCoord;","   vColor = aColor;","}"]},d.StripShader.prototype.init=function(){var a=d.gl,b=d.compileProgram(a,this.vertexSrc,this.fragmentSrc);a.useProgram(b),this.uSampler=a.getUniformLocation(b,"uSampler"),this.projectionVector=a.getUniformLocation(b,"projectionVector"),this.offsetVector=a.getUniformLocation(b,"offsetVector"),this.colorAttribute=a.getAttribLocation(b,"aColor"),this.aVertexPosition=a.getAttribLocation(b,"aVertexPosition"),this.aTextureCoord=a.getAttribLocation(b,"aTextureCoord"),this.translationMatrix=a.getUniformLocation(b,"translationMatrix"),this.alpha=a.getUniformLocation(b,"alpha"),this.program=b},d.PrimitiveShader=function(a){this.gl=a,this.program=null,this.fragmentSrc=["precision mediump float;","varying vec4 vColor;","void main(void) {","   gl_FragColor = vColor;","}"],this.vertexSrc=["attribute vec2 aVertexPosition;","attribute vec4 aColor;","uniform mat3 translationMatrix;","uniform vec2 projectionVector;","uniform vec2 offsetVector;","uniform float alpha;","uniform vec3 tint;","varying vec4 vColor;","void main(void) {","   vec3 v = translationMatrix * vec3(aVertexPosition , 1.0);","   v -= offsetVector.xyx;","   gl_Position = vec4( v.x / projectionVector.x -1.0, v.y / -projectionVector.y + 1.0 , 0.0, 1.0);","   vColor = aColor * vec4(tint * alpha, alpha);","}"],this.init()},d.PrimitiveShader.prototype.init=function(){var a=this.gl,b=d.compileProgram(a,this.vertexSrc,this.fragmentSrc);a.useProgram(b),this.projectionVector=a.getUniformLocation(b,"projectionVector"),this.offsetVector=a.getUniformLocation(b,"offsetVector"),this.tintColor=a.getUniformLocation(b,"tint"),this.aVertexPosition=a.getAttribLocation(b,"aVertexPosition"),this.colorAttribute=a.getAttribLocation(b,"aColor"),this.attributes=[this.aVertexPosition,this.colorAttribute],this.translationMatrix=a.getUniformLocation(b,"translationMatrix"),this.alpha=a.getUniformLocation(b,"alpha"),this.program=b},d.PrimitiveShader.prototype.destroy=function(){this.gl.deleteProgram(this.program),this.uniforms=null,this.gl=null,this.attribute=null},d.WebGLGraphics=function(){},d.WebGLGraphics.renderGraphics=function(a,b){var c=b.gl,e=b.projection,f=b.offset,g=b.shaderManager.primitiveShader;a._webGL[c.id]||(a._webGL[c.id]={points:[],indices:[],lastIndex:0,buffer:c.createBuffer(),indexBuffer:c.createBuffer()});var h=a._webGL[c.id];a.dirty&&(a.dirty=!1,a.clearDirty&&(a.clearDirty=!1,h.lastIndex=0,h.points=[],h.indices=[]),d.WebGLGraphics.updateGraphics(a,c)),b.shaderManager.activatePrimitiveShader(),c.blendFunc(c.ONE,c.ONE_MINUS_SRC_ALPHA),c.uniformMatrix3fv(g.translationMatrix,!1,a.worldTransform.toArray(!0)),c.uniform2f(g.projectionVector,e.x,-e.y),c.uniform2f(g.offsetVector,-f.x,-f.y),c.uniform3fv(g.tintColor,d.hex2rgb(a.tint)),c.uniform1f(g.alpha,a.worldAlpha),c.bindBuffer(c.ARRAY_BUFFER,h.buffer),c.vertexAttribPointer(g.aVertexPosition,2,c.FLOAT,!1,24,0),c.vertexAttribPointer(g.colorAttribute,4,c.FLOAT,!1,24,8),c.bindBuffer(c.ELEMENT_ARRAY_BUFFER,h.indexBuffer),c.drawElements(c.TRIANGLE_STRIP,h.indices.length,c.UNSIGNED_SHORT,0),b.shaderManager.deactivatePrimitiveShader()},d.WebGLGraphics.updateGraphics=function(a,b){for(var c=a._webGL[b.id],e=c.lastIndex;e<a.graphicsData.length;e++){var f=a.graphicsData[e];f.type===d.Graphics.POLY?(f.fill&&f.points.length>3&&d.WebGLGraphics.buildPoly(f,c),f.lineWidth>0&&d.WebGLGraphics.buildLine(f,c)):f.type===d.Graphics.RECT?d.WebGLGraphics.buildRectangle(f,c):(f.type===d.Graphics.CIRC||f.type===d.Graphics.ELIP)&&d.WebGLGraphics.buildCircle(f,c)}c.lastIndex=a.graphicsData.length,c.glPoints=new Float32Array(c.points),b.bindBuffer(b.ARRAY_BUFFER,c.buffer),b.bufferData(b.ARRAY_BUFFER,c.glPoints,b.STATIC_DRAW),c.glIndicies=new Uint16Array(c.indices),b.bindBuffer(b.ELEMENT_ARRAY_BUFFER,c.indexBuffer),b.bufferData(b.ELEMENT_ARRAY_BUFFER,c.glIndicies,b.STATIC_DRAW)},d.WebGLGraphics.buildRectangle=function(a,b){var c=a.points,e=c[0],f=c[1],g=c[2],h=c[3];if(a.fill){var i=d.hex2rgb(a.fillColor),j=a.fillAlpha,k=i[0]*j,l=i[1]*j,m=i[2]*j,n=b.points,o=b.indices,p=n.length/6;n.push(e,f),n.push(k,l,m,j),n.push(e+g,f),n.push(k,l,m,j),n.push(e,f+h),n.push(k,l,m,j),n.push(e+g,f+h),n.push(k,l,m,j),o.push(p,p,p+1,p+2,p+3,p+3)}if(a.lineWidth){var q=a.points;a.points=[e,f,e+g,f,e+g,f+h,e,f+h,e,f],d.WebGLGraphics.buildLine(a,b),a.points=q}},d.WebGLGraphics.buildCircle=function(a,b){var c=a.points,e=c[0],f=c[1],g=c[2],h=c[3],i=40,j=2*Math.PI/i,k=0;if(a.fill){var l=d.hex2rgb(a.fillColor),m=a.fillAlpha,n=l[0]*m,o=l[1]*m,p=l[2]*m,q=b.points,r=b.indices,s=q.length/6;for(r.push(s),k=0;i+1>k;k++)q.push(e,f,n,o,p,m),q.push(e+Math.sin(j*k)*g,f+Math.cos(j*k)*h,n,o,p,m),r.push(s++,s++);r.push(s-1)}if(a.lineWidth){var t=a.points;for(a.points=[],k=0;i+1>k;k++)a.points.push(e+Math.sin(j*k)*g,f+Math.cos(j*k)*h);d.WebGLGraphics.buildLine(a,b),a.points=t}},d.WebGLGraphics.buildLine=function(a,b){var c=0,e=a.points;if(0!==e.length){if(a.lineWidth%2)for(c=0;c<e.length;c++)e[c]+=.5;var f=new d.Point(e[0],e[1]),g=new d.Point(e[e.length-2],e[e.length-1]);if(f.x===g.x&&f.y===g.y){e.pop(),e.pop(),g=new d.Point(e[e.length-2],e[e.length-1]);var h=g.x+.5*(f.x-g.x),i=g.y+.5*(f.y-g.y);e.unshift(h,i),e.push(h,i)}var j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D,E,F,G=b.points,H=b.indices,I=e.length/2,J=e.length,K=G.length/6,L=a.lineWidth/2,M=d.hex2rgb(a.lineColor),N=a.lineAlpha,O=M[0]*N,P=M[1]*N,Q=M[2]*N;for(l=e[0],m=e[1],n=e[2],o=e[3],r=-(m-o),s=l-n,F=Math.sqrt(r*r+s*s),r/=F,s/=F,r*=L,s*=L,G.push(l-r,m-s,O,P,Q,N),G.push(l+r,m+s,O,P,Q,N),c=1;I-1>c;c++)l=e[2*(c-1)],m=e[2*(c-1)+1],n=e[2*c],o=e[2*c+1],p=e[2*(c+1)],q=e[2*(c+1)+1],r=-(m-o),s=l-n,F=Math.sqrt(r*r+s*s),r/=F,s/=F,r*=L,s*=L,t=-(o-q),u=n-p,F=Math.sqrt(t*t+u*u),t/=F,u/=F,t*=L,u*=L,x=-s+m-(-s+o),y=-r+n-(-r+l),z=(-r+l)*(-s+o)-(-r+n)*(-s+m),A=-u+q-(-u+o),B=-t+n-(-t+p),C=(-t+p)*(-u+o)-(-t+n)*(-u+q),D=x*B-A*y,Math.abs(D)<.1?(D+=10.1,G.push(n-r,o-s,O,P,Q,N),G.push(n+r,o+s,O,P,Q,N)):(j=(y*C-B*z)/D,k=(A*z-x*C)/D,E=(j-n)*(j-n)+(k-o)+(k-o),E>19600?(v=r-t,w=s-u,F=Math.sqrt(v*v+w*w),v/=F,w/=F,v*=L,w*=L,G.push(n-v,o-w),G.push(O,P,Q,N),G.push(n+v,o+w),G.push(O,P,Q,N),G.push(n-v,o-w),G.push(O,P,Q,N),J++):(G.push(j,k),G.push(O,P,Q,N),G.push(n-(j-n),o-(k-o)),G.push(O,P,Q,N)));for(l=e[2*(I-2)],m=e[2*(I-2)+1],n=e[2*(I-1)],o=e[2*(I-1)+1],r=-(m-o),s=l-n,F=Math.sqrt(r*r+s*s),r/=F,s/=F,r*=L,s*=L,G.push(n-r,o-s),G.push(O,P,Q,N),G.push(n+r,o+s),G.push(O,P,Q,N),H.push(K),c=0;J>c;c++)H.push(K++);H.push(K-1)}},d.WebGLGraphics.buildPoly=function(a,b){var c=a.points;if(!(c.length<6)){var e=b.points,f=b.indices,g=c.length/2,h=d.hex2rgb(a.fillColor),i=a.fillAlpha,j=h[0]*i,k=h[1]*i,l=h[2]*i,m=d.PolyK.Triangulate(c),n=e.length/6,o=0;for(o=0;o<m.length;o+=3)f.push(m[o]+n),f.push(m[o]+n),f.push(m[o+1]+n),f.push(m[o+2]+n),f.push(m[o+2]+n);for(o=0;g>o;o++)e.push(c[2*o],c[2*o+1],j,k,l,i)}},d.glContexts=[],d.WebGLRenderer=function(a,b,c,e,f){d.defaultRenderer||(d.defaultRenderer=this),this.type=d.WEBGL_RENDERER,this.transparent=!!e,this.width=a||800,this.height=b||600,this.view=c||document.createElement("canvas"),this.view.width=this.width,this.view.height=this.height,this.contextLost=this.handleContextLost.bind(this),this.contextRestoredLost=this.handleContextRestored.bind(this),this.view.addEventListener("webglcontextlost",this.contextLost,!1),this.view.addEventListener("webglcontextrestored",this.contextRestoredLost,!1),this.options={alpha:this.transparent,antialias:!!f,premultipliedAlpha:!!e,stencil:!0};try{this.gl=this.view.getContext("experimental-webgl",this.options)}catch(g){try{this.gl=this.view.getContext("webgl",this.options)}catch(h){throw new Error(" This browser does not support webGL. Try using the canvas renderer"+this)}}var i=this.gl;this.glContextId=i.id=d.WebGLRenderer.glContextId++,d.glContexts[this.glContextId]=i,d.blendModesWebGL||(d.blendModesWebGL=[],d.blendModesWebGL[d.blendModes.NORMAL]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.ADD]=[i.SRC_ALPHA,i.DST_ALPHA],d.blendModesWebGL[d.blendModes.MULTIPLY]=[i.DST_COLOR,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.SCREEN]=[i.SRC_ALPHA,i.ONE],d.blendModesWebGL[d.blendModes.OVERLAY]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.DARKEN]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.LIGHTEN]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.COLOR_DODGE]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.COLOR_BURN]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.HARD_LIGHT]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.SOFT_LIGHT]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.DIFFERENCE]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.EXCLUSION]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.HUE]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.SATURATION]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.COLOR]=[i.ONE,i.ONE_MINUS_SRC_ALPHA],d.blendModesWebGL[d.blendModes.LUMINOSITY]=[i.ONE,i.ONE_MINUS_SRC_ALPHA]),this.projection=new d.Point,this.projection.x=this.width/2,this.projection.y=-this.height/2,this.offset=new d.Point(0,0),this.resize(this.width,this.height),this.contextLost=!1,this.shaderManager=new d.WebGLShaderManager(i),this.spriteBatch=new d.WebGLSpriteBatch(i),this.maskManager=new d.WebGLMaskManager(i),this.filterManager=new d.WebGLFilterManager(i,this.transparent),this.renderSession={},this.renderSession.gl=this.gl,this.renderSession.drawCount=0,this.renderSession.shaderManager=this.shaderManager,this.renderSession.maskManager=this.maskManager,this.renderSession.filterManager=this.filterManager,this.renderSession.spriteBatch=this.spriteBatch,i.useProgram(this.shaderManager.defaultShader.program),i.disable(i.DEPTH_TEST),i.disable(i.CULL_FACE),i.enable(i.BLEND),i.colorMask(!0,!0,!0,this.transparent)},d.WebGLRenderer.prototype.constructor=d.WebGLRenderer,d.WebGLRenderer.prototype.render=function(a){if(!this.contextLost){this.__stage!==a&&(a.interactive&&a.interactionManager.removeEvents(),this.__stage=a),d.WebGLRenderer.updateTextures(),a.updateTransform();var b=this.gl;b.viewport(0,0,this.width,this.height),b.bindFramebuffer(b.FRAMEBUFFER,null),this.transparent?b.clearColor(0,0,0,0):b.clearColor(a.backgroundColorSplit[0],a.backgroundColorSplit[1],a.backgroundColorSplit[2],1),b.clear(b.COLOR_BUFFER_BIT),this.renderDisplayObject(a,this.projection),a.interactive?a._interactiveEventsAdded||(a._interactiveEventsAdded=!0,a.interactionManager.setTarget(this)):a._interactiveEventsAdded&&(a._interactiveEventsAdded=!1,a.interactionManager.setTarget(this))}},d.WebGLRenderer.prototype.renderDisplayObject=function(a,b,c){this.renderSession.drawCount=0,this.renderSession.currentBlendMode=9999,this.renderSession.projection=b,this.renderSession.offset=this.offset,this.spriteBatch.begin(this.renderSession),this.filterManager.begin(this.renderSession,c),a._renderWebGL(this.renderSession),this.spriteBatch.end()},d.WebGLRenderer.updateTextures=function(){var a=0;for(a=0;a<d.Texture.frameUpdates.length;a++)d.WebGLRenderer.updateTextureFrame(d.Texture.frameUpdates[a]);for(a=0;a<d.texturesToDestroy.length;a++)d.WebGLRenderer.destroyTexture(d.texturesToDestroy[a]);d.texturesToUpdate.length=0,d.texturesToDestroy.length=0,d.Texture.frameUpdates.length=0},d.WebGLRenderer.destroyTexture=function(a){for(var b=a._glTextures.length-1;b>=0;b--){var c=a._glTextures[b],e=d.glContexts[b];e&&c&&e.deleteTexture(c)}a._glTextures.length=0},d.WebGLRenderer.updateTextureFrame=function(a){a.updateFrame=!1,a._updateWebGLuvs()},d.WebGLRenderer.prototype.resize=function(a,b){this.width=a,this.height=b,this.view.width=a,this.view.height=b,this.gl.viewport(0,0,this.width,this.height),this.projection.x=this.width/2,this.projection.y=-this.height/2},d.createWebGLTexture=function(a,b){return a.hasLoaded&&(a._glTextures[b.id]=b.createTexture(),b.bindTexture(b.TEXTURE_2D,a._glTextures[b.id]),b.pixelStorei(b.UNPACK_PREMULTIPLY_ALPHA_WEBGL,!0),b.texImage2D(b.TEXTURE_2D,0,b.RGBA,b.RGBA,b.UNSIGNED_BYTE,a.source),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_MAG_FILTER,a.scaleMode===d.scaleModes.LINEAR?b.LINEAR:b.NEAREST),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_MIN_FILTER,a.scaleMode===d.scaleModes.LINEAR?b.LINEAR:b.NEAREST),a._powerOf2?(b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_S,b.REPEAT),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_T,b.REPEAT)):(b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_S,b.CLAMP_TO_EDGE),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_T,b.CLAMP_TO_EDGE)),b.bindTexture(b.TEXTURE_2D,null)),a._glTextures[b.id]},d.updateWebGLTexture=function(a,b){a._glTextures[b.id]&&(b.bindTexture(b.TEXTURE_2D,a._glTextures[b.id]),b.pixelStorei(b.UNPACK_PREMULTIPLY_ALPHA_WEBGL,!0),b.texImage2D(b.TEXTURE_2D,0,b.RGBA,b.RGBA,b.UNSIGNED_BYTE,a.source),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_MAG_FILTER,a.scaleMode===d.scaleModes.LINEAR?b.LINEAR:b.NEAREST),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_MIN_FILTER,a.scaleMode===d.scaleModes.LINEAR?b.LINEAR:b.NEAREST),a._powerOf2?(b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_S,b.REPEAT),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_T,b.REPEAT)):(b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_S,b.CLAMP_TO_EDGE),b.texParameteri(b.TEXTURE_2D,b.TEXTURE_WRAP_T,b.CLAMP_TO_EDGE)),b.bindTexture(b.TEXTURE_2D,null))},d.WebGLRenderer.prototype.handleContextLost=function(a){a.preventDefault(),this.contextLost=!0},d.WebGLRenderer.prototype.handleContextRestored=function(){try{this.gl=this.view.getContext("experimental-webgl",this.options)}catch(a){try{this.gl=this.view.getContext("webgl",this.options)}catch(b){throw new Error(" This browser does not support webGL. Try using the canvas renderer"+this)}}var c=this.gl;c.id=d.WebGLRenderer.glContextId++,this.shaderManager.setContext(c),this.spriteBatch.setContext(c),this.maskManager.setContext(c),this.filterManager.setContext(c),this.renderSession.gl=this.gl,c.disable(c.DEPTH_TEST),c.disable(c.CULL_FACE),c.enable(c.BLEND),c.colorMask(!0,!0,!0,this.transparent),this.gl.viewport(0,0,this.width,this.height);for(var e in d.TextureCache){var f=d.TextureCache[e].baseTexture;f._glTextures=[]}this.contextLost=!1},d.WebGLRenderer.prototype.destroy=function(){this.view.removeEventListener("webglcontextlost",this.contextLost),this.view.removeEventListener("webglcontextrestored",this.contextRestoredLost),d.glContexts[this.glContextId]=null,this.projection=null,this.offset=null,this.shaderManager.destroy(),this.spriteBatch.destroy(),this.maskManager.destroy(),this.filterManager.destroy(),this.shaderManager=null,this.spriteBatch=null,this.maskManager=null,this.filterManager=null,this.gl=null,this.renderSession=null},d.WebGLRenderer.glContextId=0,d.WebGLMaskManager=function(a){this.maskStack=[],this.maskPosition=0,this.setContext(a)},d.WebGLMaskManager.prototype.setContext=function(a){this.gl=a},d.WebGLMaskManager.prototype.pushMask=function(a,b){var c=this.gl;0===this.maskStack.length&&(c.enable(c.STENCIL_TEST),c.stencilFunc(c.ALWAYS,1,1)),this.maskStack.push(a),c.colorMask(!1,!1,!1,!0),c.stencilOp(c.KEEP,c.KEEP,c.INCR),d.WebGLGraphics.renderGraphics(a,b),c.colorMask(!0,!0,!0,!0),c.stencilFunc(c.NOTEQUAL,0,this.maskStack.length),c.stencilOp(c.KEEP,c.KEEP,c.KEEP)},d.WebGLMaskManager.prototype.popMask=function(a){var b=this.gl,c=this.maskStack.pop();c&&(b.colorMask(!1,!1,!1,!1),b.stencilOp(b.KEEP,b.KEEP,b.DECR),d.WebGLGraphics.renderGraphics(c,a),b.colorMask(!0,!0,!0,!0),b.stencilFunc(b.NOTEQUAL,0,this.maskStack.length),b.stencilOp(b.KEEP,b.KEEP,b.KEEP)),0===this.maskStack.length&&b.disable(b.STENCIL_TEST)},d.WebGLMaskManager.prototype.destroy=function(){this.maskStack=null,this.gl=null},d.WebGLShaderManager=function(a){this.maxAttibs=10,this.attribState=[],this.tempAttribState=[];for(var b=0;b<this.maxAttibs;b++)this.attribState[b]=!1;this.setContext(a)},d.WebGLShaderManager.prototype.setContext=function(a){this.gl=a,this.primitiveShader=new d.PrimitiveShader(a),this.defaultShader=new d.PixiShader(a),this.fastShader=new d.PixiFastShader(a),this.activateShader(this.defaultShader)},d.WebGLShaderManager.prototype.setAttribs=function(a){var b;for(b=0;b<this.tempAttribState.length;b++)this.tempAttribState[b]=!1;for(b=0;b<a.length;b++){var c=a[b];this.tempAttribState[c]=!0}var d=this.gl;for(b=0;b<this.attribState.length;b++)this.attribState[b]!==this.tempAttribState[b]&&(this.attribState[b]=this.tempAttribState[b],this.tempAttribState[b]?d.enableVertexAttribArray(b):d.disableVertexAttribArray(b))},d.WebGLShaderManager.prototype.activateShader=function(a){this.currentShader=a,this.gl.useProgram(a.program),this.setAttribs(a.attributes)},d.WebGLShaderManager.prototype.activatePrimitiveShader=function(){var a=this.gl;a.useProgram(this.primitiveShader.program),this.setAttribs(this.primitiveShader.attributes)},d.WebGLShaderManager.prototype.deactivatePrimitiveShader=function(){var a=this.gl;a.useProgram(this.defaultShader.program),this.setAttribs(this.defaultShader.attributes)},d.WebGLShaderManager.prototype.destroy=function(){this.attribState=null,this.tempAttribState=null,this.primitiveShader.destroy(),this.defaultShader.destroy(),this.fastShader.destroy(),this.gl=null},d.WebGLSpriteBatch=function(a){this.vertSize=6,this.size=1e4;var b=4*this.size*this.vertSize,c=6*this.size;this.vertices=new Float32Array(b),this.indices=new Uint16Array(c),this.lastIndexCount=0;for(var d=0,e=0;c>d;d+=6,e+=4)this.indices[d+0]=e+0,this.indices[d+1]=e+1,this.indices[d+2]=e+2,this.indices[d+3]=e+0,this.indices[d+4]=e+2,this.indices[d+5]=e+3;this.drawing=!1,this.currentBatchSize=0,this.currentBaseTexture=null,this.setContext(a)},d.WebGLSpriteBatch.prototype.setContext=function(a){this.gl=a,this.vertexBuffer=a.createBuffer(),this.indexBuffer=a.createBuffer(),a.bindBuffer(a.ELEMENT_ARRAY_BUFFER,this.indexBuffer),a.bufferData(a.ELEMENT_ARRAY_BUFFER,this.indices,a.STATIC_DRAW),a.bindBuffer(a.ARRAY_BUFFER,this.vertexBuffer),a.bufferData(a.ARRAY_BUFFER,this.vertices,a.DYNAMIC_DRAW),this.currentBlendMode=99999},d.WebGLSpriteBatch.prototype.begin=function(a){this.renderSession=a,this.shader=this.renderSession.shaderManager.defaultShader,this.start()},d.WebGLSpriteBatch.prototype.end=function(){this.flush()},d.WebGLSpriteBatch.prototype.render=function(a){(a.texture.baseTexture!==this.currentBaseTexture||this.currentBatchSize>=this.size)&&(this.flush(),this.currentBaseTexture=a.texture.baseTexture),a.blendMode!==this.currentBlendMode&&this.setBlendMode(a.blendMode);var b=a._uvs||a.texture._uvs;if(b){var c,d,e,f,g=a.worldAlpha,h=a.tint,i=this.vertices,j=a.texture.frame.width,k=a.texture.frame.height,l=a.anchor.x,m=a.anchor.y;if(a.texture.trim){var n=a.texture.trim;d=n.x-l*n.width,c=d+j,f=n.y-m*n.height,e=f+k}else c=j*(1-l),d=j*-l,e=k*(1-m),f=k*-m;var o=4*this.currentBatchSize*this.vertSize,p=a.worldTransform,q=p.a,r=p.c,s=p.b,t=p.d,u=p.tx,v=p.ty;i[o++]=q*d+s*f+u,i[o++]=t*f+r*d+v,i[o++]=b.x0,i[o++]=b.y0,i[o++]=g,i[o++]=h,i[o++]=q*c+s*f+u,i[o++]=t*f+r*c+v,i[o++]=b.x1,i[o++]=b.y1,i[o++]=g,i[o++]=h,i[o++]=q*c+s*e+u,i[o++]=t*e+r*c+v,i[o++]=b.x2,i[o++]=b.y2,i[o++]=g,i[o++]=h,i[o++]=q*d+s*e+u,i[o++]=t*e+r*d+v,i[o++]=b.x3,i[o++]=b.y3,i[o++]=g,i[o++]=h,this.currentBatchSize++}},d.WebGLSpriteBatch.prototype.renderTilingSprite=function(a){var b=a.tilingTexture;(b.baseTexture!==this.currentBaseTexture||this.currentBatchSize>=this.size)&&(this.flush(),this.currentBaseTexture=b.baseTexture),a.blendMode!==this.currentBlendMode&&this.setBlendMode(a.blendMode),a._uvs||(a._uvs=new d.TextureUvs);var c=a._uvs;a.tilePosition.x%=b.baseTexture.width,a.tilePosition.y%=b.baseTexture.height;var e=a.tilePosition.x/b.baseTexture.width,f=a.tilePosition.y/b.baseTexture.height,g=a.width/b.baseTexture.width/(a.tileScale.x*a.tileScaleOffset.x),h=a.height/b.baseTexture.height/(a.tileScale.y*a.tileScaleOffset.y);c.x0=0-e,c.y0=0-f,c.x1=1*g-e,c.y1=0-f,c.x2=1*g-e,c.y2=1*h-f,c.x3=0-e,c.y3=1*h-f;var i=a.worldAlpha,j=a.tint,k=this.vertices,l=a.width,m=a.height,n=a.anchor.x,o=a.anchor.y,p=l*(1-n),q=l*-n,r=m*(1-o),s=m*-o,t=4*this.currentBatchSize*this.vertSize,u=a.worldTransform,v=u.a,w=u.c,x=u.b,y=u.d,z=u.tx,A=u.ty;k[t++]=v*q+x*s+z,k[t++]=y*s+w*q+A,k[t++]=c.x0,k[t++]=c.y0,k[t++]=i,k[t++]=j,k[t++]=v*p+x*s+z,k[t++]=y*s+w*p+A,k[t++]=c.x1,k[t++]=c.y1,k[t++]=i,k[t++]=j,k[t++]=v*p+x*r+z,k[t++]=y*r+w*p+A,k[t++]=c.x2,k[t++]=c.y2,k[t++]=i,k[t++]=j,k[t++]=v*q+x*r+z,k[t++]=y*r+w*q+A,k[t++]=c.x3,k[t++]=c.y3,k[t++]=i,k[t++]=j,this.currentBatchSize++},d.WebGLSpriteBatch.prototype.flush=function(){if(0!==this.currentBatchSize){var a=this.gl;if(a.bindTexture(a.TEXTURE_2D,this.currentBaseTexture._glTextures[a.id]||d.createWebGLTexture(this.currentBaseTexture,a)),this.currentBatchSize>.5*this.size)a.bufferSubData(a.ARRAY_BUFFER,0,this.vertices);else{var b=this.vertices.subarray(0,4*this.currentBatchSize*this.vertSize);a.bufferSubData(a.ARRAY_BUFFER,0,b)}a.drawElements(a.TRIANGLES,6*this.currentBatchSize,a.UNSIGNED_SHORT,0),this.currentBatchSize=0,this.renderSession.drawCount++}},d.WebGLSpriteBatch.prototype.stop=function(){this.flush()},d.WebGLSpriteBatch.prototype.start=function(){var a=this.gl;a.activeTexture(a.TEXTURE0),a.bindBuffer(a.ARRAY_BUFFER,this.vertexBuffer),a.bindBuffer(a.ELEMENT_ARRAY_BUFFER,this.indexBuffer);var b=this.renderSession.projection;a.uniform2f(this.shader.projectionVector,b.x,b.y);var c=4*this.vertSize;a.vertexAttribPointer(this.shader.aVertexPosition,2,a.FLOAT,!1,c,0),a.vertexAttribPointer(this.shader.aTextureCoord,2,a.FLOAT,!1,c,8),a.vertexAttribPointer(this.shader.colorAttribute,2,a.FLOAT,!1,c,16),this.currentBlendMode!==d.blendModes.NORMAL&&this.setBlendMode(d.blendModes.NORMAL)},d.WebGLSpriteBatch.prototype.setBlendMode=function(a){this.flush(),this.currentBlendMode=a;var b=d.blendModesWebGL[this.currentBlendMode];
this.gl.blendFunc(b[0],b[1])},d.WebGLSpriteBatch.prototype.destroy=function(){this.vertices=null,this.indices=null,this.gl.deleteBuffer(this.vertexBuffer),this.gl.deleteBuffer(this.indexBuffer),this.currentBaseTexture=null,this.gl=null},d.WebGLFastSpriteBatch=function(a){this.vertSize=10,this.maxSize=6e3,this.size=this.maxSize;var b=4*this.size*this.vertSize,c=6*this.maxSize;this.vertices=new Float32Array(b),this.indices=new Uint16Array(c),this.vertexBuffer=null,this.indexBuffer=null,this.lastIndexCount=0;for(var d=0,e=0;c>d;d+=6,e+=4)this.indices[d+0]=e+0,this.indices[d+1]=e+1,this.indices[d+2]=e+2,this.indices[d+3]=e+0,this.indices[d+4]=e+2,this.indices[d+5]=e+3;this.drawing=!1,this.currentBatchSize=0,this.currentBaseTexture=null,this.currentBlendMode=0,this.renderSession=null,this.shader=null,this.matrix=null,this.setContext(a)},d.WebGLFastSpriteBatch.prototype.setContext=function(a){this.gl=a,this.vertexBuffer=a.createBuffer(),this.indexBuffer=a.createBuffer(),a.bindBuffer(a.ELEMENT_ARRAY_BUFFER,this.indexBuffer),a.bufferData(a.ELEMENT_ARRAY_BUFFER,this.indices,a.STATIC_DRAW),a.bindBuffer(a.ARRAY_BUFFER,this.vertexBuffer),a.bufferData(a.ARRAY_BUFFER,this.vertices,a.DYNAMIC_DRAW),this.currentBlendMode=99999},d.WebGLFastSpriteBatch.prototype.begin=function(a,b){this.renderSession=b,this.shader=this.renderSession.shaderManager.fastShader,this.matrix=a.worldTransform.toArray(!0),this.start()},d.WebGLFastSpriteBatch.prototype.end=function(){this.flush()},d.WebGLFastSpriteBatch.prototype.render=function(a){var b=a.children,c=b[0];if(c.texture._uvs){this.currentBaseTexture=c.texture.baseTexture,c.blendMode!==this.currentBlendMode&&this.setBlendMode(c.blendMode);for(var d=0,e=b.length;e>d;d++)this.renderSprite(b[d]);this.flush()}},d.WebGLFastSpriteBatch.prototype.renderSprite=function(a){if(a.texture.baseTexture===this.currentBaseTexture||(this.flush(),this.currentBaseTexture=a.texture.baseTexture,a.texture._uvs)){var b,c,d,e,f,g,h,i,j=this.vertices;if(b=a.texture._uvs,c=a.texture.frame.width,d=a.texture.frame.height,a.texture.trim){var k=a.texture.trim;f=k.x-a.anchor.x*k.width,e=f+a.texture.frame.width,h=k.y-a.anchor.y*k.height,g=h+a.texture.frame.height}else e=a.texture.frame.width*(1-a.anchor.x),f=a.texture.frame.width*-a.anchor.x,g=a.texture.frame.height*(1-a.anchor.y),h=a.texture.frame.height*-a.anchor.y;i=4*this.currentBatchSize*this.vertSize,j[i++]=f,j[i++]=h,j[i++]=a.position.x,j[i++]=a.position.y,j[i++]=a.scale.x,j[i++]=a.scale.y,j[i++]=a.rotation,j[i++]=b.x0,j[i++]=b.y1,j[i++]=a.alpha,j[i++]=e,j[i++]=h,j[i++]=a.position.x,j[i++]=a.position.y,j[i++]=a.scale.x,j[i++]=a.scale.y,j[i++]=a.rotation,j[i++]=b.x1,j[i++]=b.y1,j[i++]=a.alpha,j[i++]=e,j[i++]=g,j[i++]=a.position.x,j[i++]=a.position.y,j[i++]=a.scale.x,j[i++]=a.scale.y,j[i++]=a.rotation,j[i++]=b.x2,j[i++]=b.y2,j[i++]=a.alpha,j[i++]=f,j[i++]=g,j[i++]=a.position.x,j[i++]=a.position.y,j[i++]=a.scale.x,j[i++]=a.scale.y,j[i++]=a.rotation,j[i++]=b.x3,j[i++]=b.y3,j[i++]=a.alpha,this.currentBatchSize++,this.currentBatchSize>=this.size&&this.flush()}},d.WebGLFastSpriteBatch.prototype.flush=function(){if(0!==this.currentBatchSize){var a=this.gl;if(this.currentBaseTexture._glTextures[a.id]||d.createWebGLTexture(this.currentBaseTexture,a),a.bindTexture(a.TEXTURE_2D,this.currentBaseTexture._glTextures[a.id]),this.currentBatchSize>.5*this.size)a.bufferSubData(a.ARRAY_BUFFER,0,this.vertices);else{var b=this.vertices.subarray(0,4*this.currentBatchSize*this.vertSize);a.bufferSubData(a.ARRAY_BUFFER,0,b)}a.drawElements(a.TRIANGLES,6*this.currentBatchSize,a.UNSIGNED_SHORT,0),this.currentBatchSize=0,this.renderSession.drawCount++}},d.WebGLFastSpriteBatch.prototype.stop=function(){this.flush()},d.WebGLFastSpriteBatch.prototype.start=function(){var a=this.gl;a.activeTexture(a.TEXTURE0),a.bindBuffer(a.ARRAY_BUFFER,this.vertexBuffer),a.bindBuffer(a.ELEMENT_ARRAY_BUFFER,this.indexBuffer);var b=this.renderSession.projection;a.uniform2f(this.shader.projectionVector,b.x,b.y),a.uniformMatrix3fv(this.shader.uMatrix,!1,this.matrix);var c=4*this.vertSize;a.vertexAttribPointer(this.shader.aVertexPosition,2,a.FLOAT,!1,c,0),a.vertexAttribPointer(this.shader.aPositionCoord,2,a.FLOAT,!1,c,8),a.vertexAttribPointer(this.shader.aScale,2,a.FLOAT,!1,c,16),a.vertexAttribPointer(this.shader.aRotation,1,a.FLOAT,!1,c,24),a.vertexAttribPointer(this.shader.aTextureCoord,2,a.FLOAT,!1,c,28),a.vertexAttribPointer(this.shader.colorAttribute,1,a.FLOAT,!1,c,36),this.currentBlendMode!==d.blendModes.NORMAL&&this.setBlendMode(d.blendModes.NORMAL)},d.WebGLFastSpriteBatch.prototype.setBlendMode=function(a){this.flush(),this.currentBlendMode=a;var b=d.blendModesWebGL[this.currentBlendMode];this.gl.blendFunc(b[0],b[1])},d.WebGLFilterManager=function(a,b){this.transparent=b,this.filterStack=[],this.offsetX=0,this.offsetY=0,this.setContext(a)},d.WebGLFilterManager.prototype.setContext=function(a){this.gl=a,this.texturePool=[],this.initShaderBuffers()},d.WebGLFilterManager.prototype.begin=function(a,b){this.renderSession=a,this.defaultShader=a.shaderManager.defaultShader;var c=this.renderSession.projection;this.width=2*c.x,this.height=2*-c.y,this.buffer=b},d.WebGLFilterManager.prototype.pushFilter=function(a){var b=this.gl,c=this.renderSession.projection,e=this.renderSession.offset;this.filterStack.push(a);var f=a.filterPasses[0];this.offsetX+=a.target.filterArea.x,this.offsetY+=a.target.filterArea.y;var g=this.texturePool.pop();g?g.resize(this.width,this.height):g=new d.FilterTexture(this.gl,this.width,this.height),b.bindTexture(b.TEXTURE_2D,g.texture),a.target.filterArea=a.target.getBounds();var h=a.target.filterArea,i=f.padding;h.x-=i,h.y-=i,h.width+=2*i,h.height+=2*i,h.x<0&&(h.x=0),h.width>this.width&&(h.width=this.width),h.y<0&&(h.y=0),h.height>this.height&&(h.height=this.height),b.bindFramebuffer(b.FRAMEBUFFER,g.frameBuffer),b.viewport(0,0,h.width,h.height),c.x=h.width/2,c.y=-h.height/2,e.x=-h.x,e.y=-h.y,b.uniform2f(this.defaultShader.projectionVector,h.width/2,-h.height/2),b.uniform2f(this.defaultShader.offsetVector,-h.x,-h.y),b.colorMask(!0,!0,!0,!0),b.clearColor(0,0,0,0),b.clear(b.COLOR_BUFFER_BIT),a._glFilterTexture=g},d.WebGLFilterManager.prototype.popFilter=function(){var a=this.gl,b=this.filterStack.pop(),c=b.target.filterArea,e=b._glFilterTexture,f=this.renderSession.projection,g=this.renderSession.offset;if(b.filterPasses.length>1){a.viewport(0,0,c.width,c.height),a.bindBuffer(a.ARRAY_BUFFER,this.vertexBuffer),this.vertexArray[0]=0,this.vertexArray[1]=c.height,this.vertexArray[2]=c.width,this.vertexArray[3]=c.height,this.vertexArray[4]=0,this.vertexArray[5]=0,this.vertexArray[6]=c.width,this.vertexArray[7]=0,a.bufferSubData(a.ARRAY_BUFFER,0,this.vertexArray),a.bindBuffer(a.ARRAY_BUFFER,this.uvBuffer),this.uvArray[2]=c.width/this.width,this.uvArray[5]=c.height/this.height,this.uvArray[6]=c.width/this.width,this.uvArray[7]=c.height/this.height,a.bufferSubData(a.ARRAY_BUFFER,0,this.uvArray);var h=e,i=this.texturePool.pop();i||(i=new d.FilterTexture(this.gl,this.width,this.height)),a.bindFramebuffer(a.FRAMEBUFFER,i.frameBuffer),a.clear(a.COLOR_BUFFER_BIT),a.disable(a.BLEND);for(var j=0;j<b.filterPasses.length-1;j++){var k=b.filterPasses[j];a.bindFramebuffer(a.FRAMEBUFFER,i.frameBuffer),a.activeTexture(a.TEXTURE0),a.bindTexture(a.TEXTURE_2D,h.texture),this.applyFilterPass(k,c,c.width,c.height);var l=h;h=i,i=l}a.enable(a.BLEND),e=h,this.texturePool.push(i)}var m=b.filterPasses[b.filterPasses.length-1];this.offsetX-=c.x,this.offsetY-=c.y;var n=this.width,o=this.height,p=0,q=0,r=this.buffer;if(0===this.filterStack.length)a.colorMask(!0,!0,!0,this.transparent);else{var s=this.filterStack[this.filterStack.length-1];c=s.target.filterArea,n=c.width,o=c.height,p=c.x,q=c.y,r=s._glFilterTexture.frameBuffer}f.x=n/2,f.y=-o/2,g.x=p,g.y=q,c=b.target.filterArea;var t=c.x-p,u=c.y-q;a.bindBuffer(a.ARRAY_BUFFER,this.vertexBuffer),this.vertexArray[0]=t,this.vertexArray[1]=u+c.height,this.vertexArray[2]=t+c.width,this.vertexArray[3]=u+c.height,this.vertexArray[4]=t,this.vertexArray[5]=u,this.vertexArray[6]=t+c.width,this.vertexArray[7]=u,a.bufferSubData(a.ARRAY_BUFFER,0,this.vertexArray),a.bindBuffer(a.ARRAY_BUFFER,this.uvBuffer),this.uvArray[2]=c.width/this.width,this.uvArray[5]=c.height/this.height,this.uvArray[6]=c.width/this.width,this.uvArray[7]=c.height/this.height,a.bufferSubData(a.ARRAY_BUFFER,0,this.uvArray),a.viewport(0,0,n,o),a.bindFramebuffer(a.FRAMEBUFFER,r),a.activeTexture(a.TEXTURE0),a.bindTexture(a.TEXTURE_2D,e.texture),this.applyFilterPass(m,c,n,o),a.useProgram(this.defaultShader.program),a.uniform2f(this.defaultShader.projectionVector,n/2,-o/2),a.uniform2f(this.defaultShader.offsetVector,-p,-q),this.texturePool.push(e),b._glFilterTexture=null},d.WebGLFilterManager.prototype.applyFilterPass=function(a,b,c,e){var f=this.gl,g=a.shaders[f.id];g||(g=new d.PixiShader(f),g.fragmentSrc=a.fragmentSrc,g.uniforms=a.uniforms,g.init(),a.shaders[f.id]=g),f.useProgram(g.program),f.uniform2f(g.projectionVector,c/2,-e/2),f.uniform2f(g.offsetVector,0,0),a.uniforms.dimensions&&(a.uniforms.dimensions.value[0]=this.width,a.uniforms.dimensions.value[1]=this.height,a.uniforms.dimensions.value[2]=this.vertexArray[0],a.uniforms.dimensions.value[3]=this.vertexArray[5]),g.syncUniforms(),f.bindBuffer(f.ARRAY_BUFFER,this.vertexBuffer),f.vertexAttribPointer(g.aVertexPosition,2,f.FLOAT,!1,0,0),f.bindBuffer(f.ARRAY_BUFFER,this.uvBuffer),f.vertexAttribPointer(g.aTextureCoord,2,f.FLOAT,!1,0,0),f.bindBuffer(f.ARRAY_BUFFER,this.colorBuffer),f.vertexAttribPointer(g.colorAttribute,2,f.FLOAT,!1,0,0),f.bindBuffer(f.ELEMENT_ARRAY_BUFFER,this.indexBuffer),f.drawElements(f.TRIANGLES,6,f.UNSIGNED_SHORT,0),this.renderSession.drawCount++},d.WebGLFilterManager.prototype.initShaderBuffers=function(){var a=this.gl;this.vertexBuffer=a.createBuffer(),this.uvBuffer=a.createBuffer(),this.colorBuffer=a.createBuffer(),this.indexBuffer=a.createBuffer(),this.vertexArray=new Float32Array([0,0,1,0,0,1,1,1]),a.bindBuffer(a.ARRAY_BUFFER,this.vertexBuffer),a.bufferData(a.ARRAY_BUFFER,this.vertexArray,a.STATIC_DRAW),this.uvArray=new Float32Array([0,0,1,0,0,1,1,1]),a.bindBuffer(a.ARRAY_BUFFER,this.uvBuffer),a.bufferData(a.ARRAY_BUFFER,this.uvArray,a.STATIC_DRAW),this.colorArray=new Float32Array([1,16777215,1,16777215,1,16777215,1,16777215]),a.bindBuffer(a.ARRAY_BUFFER,this.colorBuffer),a.bufferData(a.ARRAY_BUFFER,this.colorArray,a.STATIC_DRAW),a.bindBuffer(a.ELEMENT_ARRAY_BUFFER,this.indexBuffer),a.bufferData(a.ELEMENT_ARRAY_BUFFER,new Uint16Array([0,1,2,1,3,2]),a.STATIC_DRAW)},d.WebGLFilterManager.prototype.destroy=function(){var a=this.gl;this.filterStack=null,this.offsetX=0,this.offsetY=0;for(var b=0;b<this.texturePool.length;b++)this.texturePool.destroy();this.texturePool=null,a.deleteBuffer(this.vertexBuffer),a.deleteBuffer(this.uvBuffer),a.deleteBuffer(this.colorBuffer),a.deleteBuffer(this.indexBuffer)},d.FilterTexture=function(a,b,c){this.gl=a,this.frameBuffer=a.createFramebuffer(),this.texture=a.createTexture(),a.bindTexture(a.TEXTURE_2D,this.texture),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MAG_FILTER,a.LINEAR),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MIN_FILTER,a.LINEAR),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_S,a.CLAMP_TO_EDGE),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_T,a.CLAMP_TO_EDGE),a.bindFramebuffer(a.FRAMEBUFFER,this.framebuffer),a.bindFramebuffer(a.FRAMEBUFFER,this.frameBuffer),a.framebufferTexture2D(a.FRAMEBUFFER,a.COLOR_ATTACHMENT0,a.TEXTURE_2D,this.texture,0),this.resize(b,c)},d.FilterTexture.prototype.clear=function(){var a=this.gl;a.clearColor(0,0,0,0),a.clear(a.COLOR_BUFFER_BIT)},d.FilterTexture.prototype.resize=function(a,b){if(this.width!==a||this.height!==b){this.width=a,this.height=b;var c=this.gl;c.bindTexture(c.TEXTURE_2D,this.texture),c.texImage2D(c.TEXTURE_2D,0,c.RGBA,a,b,0,c.RGBA,c.UNSIGNED_BYTE,null)}},d.FilterTexture.prototype.destroy=function(){var a=this.gl;a.deleteFramebuffer(this.frameBuffer),a.deleteTexture(this.texture),this.frameBuffer=null,this.texture=null},d.CanvasMaskManager=function(){},d.CanvasMaskManager.prototype.pushMask=function(a,b){b.save();var c=a.alpha,e=a.worldTransform;b.setTransform(e.a,e.c,e.b,e.d,e.tx,e.ty),d.CanvasGraphics.renderGraphicsMask(a,b),b.clip(),a.worldAlpha=c},d.CanvasMaskManager.prototype.popMask=function(a){a.restore()},d.CanvasTinter=function(){},d.CanvasTinter.getTintedTexture=function(a,b){var c=a.texture;b=d.CanvasTinter.roundColor(b);var e="#"+("00000"+(0|b).toString(16)).substr(-6);if(c.tintCache=c.tintCache||{},c.tintCache[e])return c.tintCache[e];var f=d.CanvasTinter.canvas||document.createElement("canvas");if(d.CanvasTinter.tintMethod(c,b,f),d.CanvasTinter.convertTintToImage){var g=new Image;g.src=f.toDataURL(),c.tintCache[e]=g}else c.tintCache[e]=f,d.CanvasTinter.canvas=null;return f},d.CanvasTinter.tintWithMultiply=function(a,b,c){var d=c.getContext("2d"),e=a.frame;c.width=e.width,c.height=e.height,d.fillStyle="#"+("00000"+(0|b).toString(16)).substr(-6),d.fillRect(0,0,e.width,e.height),d.globalCompositeOperation="multiply",d.drawImage(a.baseTexture.source,e.x,e.y,e.width,e.height,0,0,e.width,e.height),d.globalCompositeOperation="destination-atop",d.drawImage(a.baseTexture.source,e.x,e.y,e.width,e.height,0,0,e.width,e.height)},d.CanvasTinter.tintWithOverlay=function(a,b,c){var d=c.getContext("2d"),e=a.frame;c.width=e.width,c.height=e.height,d.globalCompositeOperation="copy",d.fillStyle="#"+("00000"+(0|b).toString(16)).substr(-6),d.fillRect(0,0,e.width,e.height),d.globalCompositeOperation="destination-atop",d.drawImage(a.baseTexture.source,e.x,e.y,e.width,e.height,0,0,e.width,e.height)},d.CanvasTinter.tintWithPerPixel=function(a,b,c){var e=c.getContext("2d"),f=a.frame;c.width=f.width,c.height=f.height,e.globalCompositeOperation="copy",e.drawImage(a.baseTexture.source,f.x,f.y,f.width,f.height,0,0,f.width,f.height);for(var g=d.hex2rgb(b),h=g[0],i=g[1],j=g[2],k=e.getImageData(0,0,f.width,f.height),l=k.data,m=0;m<l.length;m+=4)l[m+0]*=h,l[m+1]*=i,l[m+2]*=j;e.putImageData(k,0,0)},d.CanvasTinter.roundColor=function(a){var b=d.CanvasTinter.cacheStepsPerColorChannel,c=d.hex2rgb(a);return c[0]=Math.min(255,c[0]/b*b),c[1]=Math.min(255,c[1]/b*b),c[2]=Math.min(255,c[2]/b*b),d.rgb2hex(c)},d.CanvasTinter.cacheStepsPerColorChannel=8,d.CanvasTinter.convertTintToImage=!1,d.CanvasTinter.canUseMultiply=d.canUseNewCanvasBlendModes(),d.CanvasTinter.tintMethod=d.CanvasTinter.canUseMultiply?d.CanvasTinter.tintWithMultiply:d.CanvasTinter.tintWithPerPixel,d.CanvasRenderer=function(a,b,c,e){d.defaultRenderer=d.defaultRenderer||this,this.type=d.CANVAS_RENDERER,this.clearBeforeRender=!0,this.roundPixels=!1,this.transparent=!!e,d.blendModesCanvas||(d.blendModesCanvas=[],d.canUseNewCanvasBlendModes()?(d.blendModesCanvas[d.blendModes.NORMAL]="source-over",d.blendModesCanvas[d.blendModes.ADD]="lighter",d.blendModesCanvas[d.blendModes.MULTIPLY]="multiply",d.blendModesCanvas[d.blendModes.SCREEN]="screen",d.blendModesCanvas[d.blendModes.OVERLAY]="overlay",d.blendModesCanvas[d.blendModes.DARKEN]="darken",d.blendModesCanvas[d.blendModes.LIGHTEN]="lighten",d.blendModesCanvas[d.blendModes.COLOR_DODGE]="color-dodge",d.blendModesCanvas[d.blendModes.COLOR_BURN]="color-burn",d.blendModesCanvas[d.blendModes.HARD_LIGHT]="hard-light",d.blendModesCanvas[d.blendModes.SOFT_LIGHT]="soft-light",d.blendModesCanvas[d.blendModes.DIFFERENCE]="difference",d.blendModesCanvas[d.blendModes.EXCLUSION]="exclusion",d.blendModesCanvas[d.blendModes.HUE]="hue",d.blendModesCanvas[d.blendModes.SATURATION]="saturation",d.blendModesCanvas[d.blendModes.COLOR]="color",d.blendModesCanvas[d.blendModes.LUMINOSITY]="luminosity"):(d.blendModesCanvas[d.blendModes.NORMAL]="source-over",d.blendModesCanvas[d.blendModes.ADD]="lighter",d.blendModesCanvas[d.blendModes.MULTIPLY]="source-over",d.blendModesCanvas[d.blendModes.SCREEN]="source-over",d.blendModesCanvas[d.blendModes.OVERLAY]="source-over",d.blendModesCanvas[d.blendModes.DARKEN]="source-over",d.blendModesCanvas[d.blendModes.LIGHTEN]="source-over",d.blendModesCanvas[d.blendModes.COLOR_DODGE]="source-over",d.blendModesCanvas[d.blendModes.COLOR_BURN]="source-over",d.blendModesCanvas[d.blendModes.HARD_LIGHT]="source-over",d.blendModesCanvas[d.blendModes.SOFT_LIGHT]="source-over",d.blendModesCanvas[d.blendModes.DIFFERENCE]="source-over",d.blendModesCanvas[d.blendModes.EXCLUSION]="source-over",d.blendModesCanvas[d.blendModes.HUE]="source-over",d.blendModesCanvas[d.blendModes.SATURATION]="source-over",d.blendModesCanvas[d.blendModes.COLOR]="source-over",d.blendModesCanvas[d.blendModes.LUMINOSITY]="source-over")),this.width=a||800,this.height=b||600,this.view=c||document.createElement("canvas"),this.context=this.view.getContext("2d",{alpha:this.transparent}),this.refresh=!0,this.view.width=this.width,this.view.height=this.height,this.count=0,this.maskManager=new d.CanvasMaskManager,this.renderSession={context:this.context,maskManager:this.maskManager,scaleMode:null,smoothProperty:null},"imageSmoothingEnabled"in this.context?this.renderSession.smoothProperty="imageSmoothingEnabled":"webkitImageSmoothingEnabled"in this.context?this.renderSession.smoothProperty="webkitImageSmoothingEnabled":"mozImageSmoothingEnabled"in this.context?this.renderSession.smoothProperty="mozImageSmoothingEnabled":"oImageSmoothingEnabled"in this.context&&(this.renderSession.smoothProperty="oImageSmoothingEnabled")},d.CanvasRenderer.prototype.constructor=d.CanvasRenderer,d.CanvasRenderer.prototype.render=function(a){d.texturesToUpdate.length=0,d.texturesToDestroy.length=0,a.updateTransform(),this.context.setTransform(1,0,0,1,0,0),this.context.globalAlpha=1,!this.transparent&&this.clearBeforeRender?(this.context.fillStyle=a.backgroundColorString,this.context.fillRect(0,0,this.width,this.height)):this.transparent&&this.clearBeforeRender&&this.context.clearRect(0,0,this.width,this.height),this.renderDisplayObject(a),a.interactive&&(a._interactiveEventsAdded||(a._interactiveEventsAdded=!0,a.interactionManager.setTarget(this))),d.Texture.frameUpdates.length>0&&(d.Texture.frameUpdates.length=0)},d.CanvasRenderer.prototype.resize=function(a,b){this.width=a,this.height=b,this.view.width=a,this.view.height=b},d.CanvasRenderer.prototype.renderDisplayObject=function(a,b){this.renderSession.context=b||this.context,a._renderCanvas(this.renderSession)},d.CanvasRenderer.prototype.renderStripFlat=function(a){var b=this.context,c=a.verticies,d=c.length/2;this.count++,b.beginPath();for(var e=1;d-2>e;e++){var f=2*e,g=c[f],h=c[f+2],i=c[f+4],j=c[f+1],k=c[f+3],l=c[f+5];b.moveTo(g,j),b.lineTo(h,k),b.lineTo(i,l)}b.fillStyle="#FF0000",b.fill(),b.closePath()},d.CanvasRenderer.prototype.renderStrip=function(a){var b=this.context,c=a.verticies,d=a.uvs,e=c.length/2;this.count++;for(var f=1;e-2>f;f++){var g=2*f,h=c[g],i=c[g+2],j=c[g+4],k=c[g+1],l=c[g+3],m=c[g+5],n=d[g]*a.texture.width,o=d[g+2]*a.texture.width,p=d[g+4]*a.texture.width,q=d[g+1]*a.texture.height,r=d[g+3]*a.texture.height,s=d[g+5]*a.texture.height;b.save(),b.beginPath(),b.moveTo(h,k),b.lineTo(i,l),b.lineTo(j,m),b.closePath(),b.clip();var t=n*r+q*p+o*s-r*p-q*o-n*s,u=h*r+q*j+i*s-r*j-q*i-h*s,v=n*i+h*p+o*j-i*p-h*o-n*j,w=n*r*j+q*i*p+h*o*s-h*r*p-q*o*j-n*i*s,x=k*r+q*m+l*s-r*m-q*l-k*s,y=n*l+k*p+o*m-l*p-k*o-n*m,z=n*r*m+q*l*p+k*o*s-k*r*p-q*o*m-n*l*s;b.transform(u/t,x/t,v/t,y/t,w/t,z/t),b.drawImage(a.texture.baseTexture.source,0,0),b.restore()}},d.CanvasBuffer=function(a,b){this.width=a,this.height=b,this.canvas=document.createElement("canvas"),this.context=this.canvas.getContext("2d"),this.canvas.width=a,this.canvas.height=b},d.CanvasBuffer.prototype.clear=function(){this.context.clearRect(0,0,this.width,this.height)},d.CanvasBuffer.prototype.resize=function(a,b){this.width=this.canvas.width=a,this.height=this.canvas.height=b},d.CanvasGraphics=function(){},d.CanvasGraphics.renderGraphics=function(a,b){for(var c=a.worldAlpha,e="",f=0;f<a.graphicsData.length;f++){var g=a.graphicsData[f],h=g.points;if(b.strokeStyle=e="#"+("00000"+(0|g.lineColor).toString(16)).substr(-6),b.lineWidth=g.lineWidth,g.type===d.Graphics.POLY){b.beginPath(),b.moveTo(h[0],h[1]);for(var i=1;i<h.length/2;i++)b.lineTo(h[2*i],h[2*i+1]);h[0]===h[h.length-2]&&h[1]===h[h.length-1]&&b.closePath(),g.fill&&(b.globalAlpha=g.fillAlpha*c,b.fillStyle=e="#"+("00000"+(0|g.fillColor).toString(16)).substr(-6),b.fill()),g.lineWidth&&(b.globalAlpha=g.lineAlpha*c,b.stroke())}else if(g.type===d.Graphics.RECT)(g.fillColor||0===g.fillColor)&&(b.globalAlpha=g.fillAlpha*c,b.fillStyle=e="#"+("00000"+(0|g.fillColor).toString(16)).substr(-6),b.fillRect(h[0],h[1],h[2],h[3])),g.lineWidth&&(b.globalAlpha=g.lineAlpha*c,b.strokeRect(h[0],h[1],h[2],h[3]));else if(g.type===d.Graphics.CIRC)b.beginPath(),b.arc(h[0],h[1],h[2],0,2*Math.PI),b.closePath(),g.fill&&(b.globalAlpha=g.fillAlpha*c,b.fillStyle=e="#"+("00000"+(0|g.fillColor).toString(16)).substr(-6),b.fill()),g.lineWidth&&(b.globalAlpha=g.lineAlpha*c,b.stroke());else if(g.type===d.Graphics.ELIP){var j=g.points,k=2*j[2],l=2*j[3],m=j[0]-k/2,n=j[1]-l/2;b.beginPath();var o=.5522848,p=k/2*o,q=l/2*o,r=m+k,s=n+l,t=m+k/2,u=n+l/2;b.moveTo(m,u),b.bezierCurveTo(m,u-q,t-p,n,t,n),b.bezierCurveTo(t+p,n,r,u-q,r,u),b.bezierCurveTo(r,u+q,t+p,s,t,s),b.bezierCurveTo(t-p,s,m,u+q,m,u),b.closePath(),g.fill&&(b.globalAlpha=g.fillAlpha*c,b.fillStyle=e="#"+("00000"+(0|g.fillColor).toString(16)).substr(-6),b.fill()),g.lineWidth&&(b.globalAlpha=g.lineAlpha*c,b.stroke())}}},d.CanvasGraphics.renderGraphicsMask=function(a,b){var c=a.graphicsData.length;if(0!==c){c>1&&(c=1,window.console.log("Pixi.js warning: masks in canvas can only mask using the first path in the graphics object"));for(var e=0;1>e;e++){var f=a.graphicsData[e],g=f.points;if(f.type===d.Graphics.POLY){b.beginPath(),b.moveTo(g[0],g[1]);for(var h=1;h<g.length/2;h++)b.lineTo(g[2*h],g[2*h+1]);g[0]===g[g.length-2]&&g[1]===g[g.length-1]&&b.closePath()}else if(f.type===d.Graphics.RECT)b.beginPath(),b.rect(g[0],g[1],g[2],g[3]),b.closePath();else if(f.type===d.Graphics.CIRC)b.beginPath(),b.arc(g[0],g[1],g[2],0,2*Math.PI),b.closePath();else if(f.type===d.Graphics.ELIP){var i=f.points,j=2*i[2],k=2*i[3],l=i[0]-j/2,m=i[1]-k/2;b.beginPath();var n=.5522848,o=j/2*n,p=k/2*n,q=l+j,r=m+k,s=l+j/2,t=m+k/2;b.moveTo(l,t),b.bezierCurveTo(l,t-p,s-o,m,s,m),b.bezierCurveTo(s+o,m,q,t-p,q,t),b.bezierCurveTo(q,t+p,s+o,r,s,r),b.bezierCurveTo(s-o,r,l,t+p,l,t),b.closePath()}}}},d.Graphics=function(){d.DisplayObjectContainer.call(this),this.renderable=!0,this.fillAlpha=1,this.lineWidth=0,this.lineColor="black",this.graphicsData=[],this.tint=16777215,this.blendMode=d.blendModes.NORMAL,this.currentPath={points:[]},this._webGL=[],this.isMask=!1,this.bounds=null,this.boundsPadding=10},d.Graphics.prototype=Object.create(d.DisplayObjectContainer.prototype),d.Graphics.prototype.constructor=d.Graphics,Object.defineProperty(d.Graphics.prototype,"cacheAsBitmap",{get:function(){return this._cacheAsBitmap},set:function(a){this._cacheAsBitmap=a,this._cacheAsBitmap?this._generateCachedSprite():(this.destroyCachedSprite(),this.dirty=!0)}}),d.Graphics.prototype.lineStyle=function(a,b,c){return this.currentPath.points.length||this.graphicsData.pop(),this.lineWidth=a||0,this.lineColor=b||0,this.lineAlpha=arguments.length<3?1:c,this.currentPath={lineWidth:this.lineWidth,lineColor:this.lineColor,lineAlpha:this.lineAlpha,fillColor:this.fillColor,fillAlpha:this.fillAlpha,fill:this.filling,points:[],type:d.Graphics.POLY},this.graphicsData.push(this.currentPath),this},d.Graphics.prototype.moveTo=function(a,b){return this.currentPath.points.length||this.graphicsData.pop(),this.currentPath=this.currentPath={lineWidth:this.lineWidth,lineColor:this.lineColor,lineAlpha:this.lineAlpha,fillColor:this.fillColor,fillAlpha:this.fillAlpha,fill:this.filling,points:[],type:d.Graphics.POLY},this.currentPath.points.push(a,b),this.graphicsData.push(this.currentPath),this},d.Graphics.prototype.lineTo=function(a,b){return this.currentPath.points.push(a,b),this.dirty=!0,this},d.Graphics.prototype.beginFill=function(a,b){return this.filling=!0,this.fillColor=a||0,this.fillAlpha=arguments.length<2?1:b,this},d.Graphics.prototype.endFill=function(){return this.filling=!1,this.fillColor=null,this.fillAlpha=1,this},d.Graphics.prototype.drawRect=function(a,b,c,e){return this.currentPath.points.length||this.graphicsData.pop(),this.currentPath={lineWidth:this.lineWidth,lineColor:this.lineColor,lineAlpha:this.lineAlpha,fillColor:this.fillColor,fillAlpha:this.fillAlpha,fill:this.filling,points:[a,b,c,e],type:d.Graphics.RECT},this.graphicsData.push(this.currentPath),this.dirty=!0,this},d.Graphics.prototype.drawCircle=function(a,b,c){return this.currentPath.points.length||this.graphicsData.pop(),this.currentPath={lineWidth:this.lineWidth,lineColor:this.lineColor,lineAlpha:this.lineAlpha,fillColor:this.fillColor,fillAlpha:this.fillAlpha,fill:this.filling,points:[a,b,c,c],type:d.Graphics.CIRC},this.graphicsData.push(this.currentPath),this.dirty=!0,this},d.Graphics.prototype.drawEllipse=function(a,b,c,e){return this.currentPath.points.length||this.graphicsData.pop(),this.currentPath={lineWidth:this.lineWidth,lineColor:this.lineColor,lineAlpha:this.lineAlpha,fillColor:this.fillColor,fillAlpha:this.fillAlpha,fill:this.filling,points:[a,b,c,e],type:d.Graphics.ELIP},this.graphicsData.push(this.currentPath),this.dirty=!0,this},d.Graphics.prototype.clear=function(){return this.lineWidth=0,this.filling=!1,this.dirty=!0,this.clearDirty=!0,this.graphicsData=[],this.bounds=null,this},d.Graphics.prototype.generateTexture=function(){var a=this.getBounds(),b=new d.CanvasBuffer(a.width,a.height),c=d.Texture.fromCanvas(b.canvas);return b.context.translate(-a.x,-a.y),d.CanvasGraphics.renderGraphics(this,b.context),c},d.Graphics.prototype._renderWebGL=function(a){if(this.visible!==!1&&0!==this.alpha&&this.isMask!==!0){if(this._cacheAsBitmap)return this.dirty&&(this._generateCachedSprite(),d.updateWebGLTexture(this._cachedSprite.texture.baseTexture,a.gl),this.dirty=!1),d.Sprite.prototype._renderWebGL.call(this._cachedSprite,a),void 0;if(a.spriteBatch.stop(),this._mask&&a.maskManager.pushMask(this.mask,a),this._filters&&a.filterManager.pushFilter(this._filterBlock),this.blendMode!==a.spriteBatch.currentBlendMode){a.spriteBatch.currentBlendMode=this.blendMode;var b=d.blendModesWebGL[a.spriteBatch.currentBlendMode];a.spriteBatch.gl.blendFunc(b[0],b[1])}if(d.WebGLGraphics.renderGraphics(this,a),this.children.length){a.spriteBatch.start();for(var c=0,e=this.children.length;e>c;c++)this.children[c]._renderWebGL(a);a.spriteBatch.stop()}this._filters&&a.filterManager.popFilter(),this._mask&&a.maskManager.popMask(a),a.drawCount++,a.spriteBatch.start()}},d.Graphics.prototype._renderCanvas=function(a){if(this.visible!==!1&&0!==this.alpha&&this.isMask!==!0){var b=a.context,c=this.worldTransform;this.blendMode!==a.currentBlendMode&&(a.currentBlendMode=this.blendMode,b.globalCompositeOperation=d.blendModesCanvas[a.currentBlendMode]),b.setTransform(c.a,c.c,c.b,c.d,c.tx,c.ty),d.CanvasGraphics.renderGraphics(this,b);for(var e=0,f=this.children.length;f>e;e++)this.children[e]._renderCanvas(a)}},d.Graphics.prototype.getBounds=function(a){this.bounds||this.updateBounds();var b=this.bounds.x,c=this.bounds.width+this.bounds.x,d=this.bounds.y,e=this.bounds.height+this.bounds.y,f=a||this.worldTransform,g=f.a,h=f.c,i=f.b,j=f.d,k=f.tx,l=f.ty,m=g*c+i*e+k,n=j*e+h*c+l,o=g*b+i*e+k,p=j*e+h*b+l,q=g*b+i*d+k,r=j*d+h*b+l,s=g*c+i*d+k,t=j*d+h*c+l,u=-1/0,v=-1/0,w=1/0,x=1/0;w=w>m?m:w,w=w>o?o:w,w=w>q?q:w,w=w>s?s:w,x=x>n?n:x,x=x>p?p:x,x=x>r?r:x,x=x>t?t:x,u=m>u?m:u,u=o>u?o:u,u=q>u?q:u,u=s>u?s:u,v=n>v?n:v,v=p>v?p:v,v=r>v?r:v,v=t>v?t:v;var y=this._bounds;return y.x=w,y.width=u-w,y.y=x,y.height=v-x,y},d.Graphics.prototype.updateBounds=function(){for(var a,b,c,e,f,g=1/0,h=-1/0,i=1/0,j=-1/0,k=0;k<this.graphicsData.length;k++){var l=this.graphicsData[k],m=l.type,n=l.lineWidth;if(a=l.points,m===d.Graphics.RECT)b=a[0]-n/2,c=a[1]-n/2,e=a[2]+n,f=a[3]+n,g=g>b?b:g,h=b+e>h?b+e:h,i=i>c?b:i,j=c+f>j?c+f:j;else if(m===d.Graphics.CIRC||m===d.Graphics.ELIP)b=a[0],c=a[1],e=a[2]+n/2,f=a[3]+n/2,g=g>b-e?b-e:g,h=b+e>h?b+e:h,i=i>c-f?c-f:i,j=c+f>j?c+f:j;else for(var o=0;o<a.length;o+=2)b=a[o],c=a[o+1],g=g>b-n?b-n:g,h=b+n>h?b+n:h,i=i>c-n?c-n:i,j=c+n>j?c+n:j}var p=this.boundsPadding;this.bounds=new d.Rectangle(g-p,i-p,h-g+2*p,j-i+2*p)},d.Graphics.prototype._generateCachedSprite=function(){var a=this.getLocalBounds();if(this._cachedSprite)this._cachedSprite.buffer.resize(a.width,a.height);else{var b=new d.CanvasBuffer(a.width,a.height),c=d.Texture.fromCanvas(b.canvas);this._cachedSprite=new d.Sprite(c),this._cachedSprite.buffer=b,this._cachedSprite.worldTransform=this.worldTransform}this._cachedSprite.anchor.x=-(a.x/a.width),this._cachedSprite.anchor.y=-(a.y/a.height),this._cachedSprite.buffer.context.translate(-a.x,-a.y),d.CanvasGraphics.renderGraphics(this,this._cachedSprite.buffer.context)},d.Graphics.prototype.destroyCachedSprite=function(){this._cachedSprite.texture.destroy(!0),this._cachedSprite=null},d.Graphics.POLY=0,d.Graphics.RECT=1,d.Graphics.CIRC=2,d.Graphics.ELIP=3,d.Strip=function(a,b,c){d.DisplayObjectContainer.call(this),this.texture=a,this.blendMode=d.blendModes.NORMAL;try{this.uvs=new Float32Array([0,1,1,1,1,0,0,1]),this.verticies=new Float32Array([0,0,0,0,0,0,0,0,0]),this.colors=new Float32Array([1,1,1,1]),this.indices=new Uint16Array([0,1,2,3])}catch(e){this.uvs=[0,1,1,1,1,0,0,1],this.verticies=[0,0,0,0,0,0,0,0,0],this.colors=[1,1,1,1],this.indices=[0,1,2,3]}this.width=b,this.height=c,a.baseTexture.hasLoaded?(this.width=this.texture.frame.width,this.height=this.texture.frame.height,this.updateFrame=!0):(this.onTextureUpdateBind=this.onTextureUpdate.bind(this),this.texture.addEventListener("update",this.onTextureUpdateBind)),this.renderable=!0},d.Strip.prototype=Object.create(d.DisplayObjectContainer.prototype),d.Strip.prototype.constructor=d.Strip,d.Strip.prototype.setTexture=function(a){this.texture=a,this.width=a.frame.width,this.height=a.frame.height,this.updateFrame=!0},d.Strip.prototype.onTextureUpdate=function(){this.updateFrame=!0},d.Rope=function(a,b){d.Strip.call(this,a),this.points=b;try{this.verticies=new Float32Array(4*b.length),this.uvs=new Float32Array(4*b.length),this.colors=new Float32Array(2*b.length),this.indices=new Uint16Array(2*b.length)}catch(c){this.verticies=new Array(4*b.length),this.uvs=new Array(4*b.length),this.colors=new Array(2*b.length),this.indices=new Array(2*b.length)}this.refresh()},d.Rope.prototype=Object.create(d.Strip.prototype),d.Rope.prototype.constructor=d.Rope,d.Rope.prototype.refresh=function(){var a=this.points;if(!(a.length<1)){var b=this.uvs,c=a[0],d=this.indices,e=this.colors;this.count-=.2,b[0]=0,b[1]=1,b[2]=0,b[3]=1,e[0]=1,e[1]=1,d[0]=0,d[1]=1;for(var f,g,h,i=a.length,j=1;i>j;j++)f=a[j],g=4*j,h=j/(i-1),j%2?(b[g]=h,b[g+1]=0,b[g+2]=h,b[g+3]=1):(b[g]=h,b[g+1]=0,b[g+2]=h,b[g+3]=1),g=2*j,e[g]=1,e[g+1]=1,g=2*j,d[g]=g,d[g+1]=g+1,c=f}},d.Rope.prototype.updateTransform=function(){var a=this.points;if(!(a.length<1)){var b,c=a[0],e={x:0,y:0};this.count-=.2;var f=this.verticies;f[0]=c.x+e.x,f[1]=c.y+e.y,f[2]=c.x-e.x,f[3]=c.y-e.y;for(var g,h,i,j,k,l=a.length,m=1;l>m;m++)g=a[m],h=4*m,b=m<a.length-1?a[m+1]:g,e.y=-(b.x-c.x),e.x=b.y-c.y,i=10*(1-m/(l-1)),i>1&&(i=1),j=Math.sqrt(e.x*e.x+e.y*e.y),k=this.texture.height/2,e.x/=j,e.y/=j,e.x*=k,e.y*=k,f[h]=g.x+e.x,f[h+1]=g.y+e.y,f[h+2]=g.x-e.x,f[h+3]=g.y-e.y,c=g;d.DisplayObjectContainer.prototype.updateTransform.call(this)}},d.Rope.prototype.setTexture=function(a){this.texture=a,this.updateFrame=!0},d.TilingSprite=function(a,b,c){d.Sprite.call(this,a),this.width=b||100,this.height=c||100,this.tileScale=new d.Point(1,1),this.tileScaleOffset=new d.Point(1,1),this.tilePosition=new d.Point(0,0),this.renderable=!0,this.tint=16777215,this.blendMode=d.blendModes.NORMAL},d.TilingSprite.prototype=Object.create(d.Sprite.prototype),d.TilingSprite.prototype.constructor=d.TilingSprite,Object.defineProperty(d.TilingSprite.prototype,"width",{get:function(){return this._width},set:function(a){this._width=a}}),Object.defineProperty(d.TilingSprite.prototype,"height",{get:function(){return this._height
},set:function(a){this._height=a}}),d.TilingSprite.prototype.onTextureUpdate=function(){this.updateFrame=!0},d.TilingSprite.prototype._renderWebGL=function(a){if(this.visible!==!1&&0!==this.alpha){var b,c;if(this.mask||this.filters){for(this.mask&&(a.spriteBatch.stop(),a.maskManager.pushMask(this.mask,a),a.spriteBatch.start()),this.filters&&(a.spriteBatch.flush(),a.filterManager.pushFilter(this._filterBlock)),this.tilingTexture?a.spriteBatch.renderTilingSprite(this):this.generateTilingTexture(!0),b=0,c=this.children.length;c>b;b++)this.children[b]._renderWebGL(a);a.spriteBatch.stop(),this.filters&&a.filterManager.popFilter(),this.mask&&a.maskManager.popMask(a),a.spriteBatch.start()}else for(this.tilingTexture?a.spriteBatch.renderTilingSprite(this):this.generateTilingTexture(!0),b=0,c=this.children.length;c>b;b++)this.children[b]._renderWebGL(a)}},d.TilingSprite.prototype._renderCanvas=function(a){if(this.visible!==!1&&0!==this.alpha){var b=a.context;this._mask&&a.maskManager.pushMask(this._mask,b),b.globalAlpha=this.worldAlpha;var c=this.worldTransform;b.setTransform(c.a,c.c,c.b,c.d,c.tx,c.ty),this.__tilePattern||(this.generateTilingTexture(!1),this.tilingTexture&&(this.__tilePattern=b.createPattern(this.tilingTexture.baseTexture.source,"repeat"))),this.blendMode!==a.currentBlendMode&&(a.currentBlendMode=this.blendMode,b.globalCompositeOperation=d.blendModesCanvas[a.currentBlendMode]),b.beginPath();var e=this.tilePosition,f=this.tileScale;e.x%=this.tilingTexture.baseTexture.width,e.y%=this.tilingTexture.baseTexture.height,b.scale(f.x,f.y),b.translate(e.x,e.y),b.fillStyle=this.__tilePattern,b.fillRect(-e.x,-e.y,this.width/f.x,this.height/f.y),b.scale(1/f.x,1/f.y),b.translate(-e.x,-e.y),b.closePath(),this._mask&&a.maskManager.popMask(a.context)}},d.TilingSprite.prototype.getBounds=function(){var a=this._width,b=this._height,c=a*(1-this.anchor.x),d=a*-this.anchor.x,e=b*(1-this.anchor.y),f=b*-this.anchor.y,g=this.worldTransform,h=g.a,i=g.c,j=g.b,k=g.d,l=g.tx,m=g.ty,n=h*d+j*f+l,o=k*f+i*d+m,p=h*c+j*f+l,q=k*f+i*c+m,r=h*c+j*e+l,s=k*e+i*c+m,t=h*d+j*e+l,u=k*e+i*d+m,v=-1/0,w=-1/0,x=1/0,y=1/0;x=x>n?n:x,x=x>p?p:x,x=x>r?r:x,x=x>t?t:x,y=y>o?o:y,y=y>q?q:y,y=y>s?s:y,y=y>u?u:y,v=n>v?n:v,v=p>v?p:v,v=r>v?r:v,v=t>v?t:v,w=o>w?o:w,w=q>w?q:w,w=s>w?s:w,w=u>w?u:w;var z=this._bounds;return z.x=x,z.width=v-x,z.y=y,z.height=w-y,this._currentBounds=z,z},d.TilingSprite.prototype.generateTilingTexture=function(a){var b=this.texture;if(b.baseTexture.hasLoaded){var c,e,f=b.baseTexture,g=b.frame,h=g.width!==f.width||g.height!==f.height;this.tilingTexture=b;var i=!1;if(a?(c=d.getNextPowerOfTwo(b.frame.width),e=d.getNextPowerOfTwo(b.frame.height),g.width!==c&&g.height!==e&&(i=!0)):h&&(c=g.width,e=g.height,i=!0),i){var j=new d.CanvasBuffer(c,e);j.context.drawImage(b.baseTexture.source,g.x,g.y,g.width,g.height,0,0,c,e),this.tilingTexture=d.Texture.fromCanvas(j.canvas),this.tileScaleOffset.x=g.width/c,this.tileScaleOffset.y=g.height/e}this.tilingTexture.baseTexture._powerOf2=!0}};var i={};i.BoneData=function(a,b){this.name=a,this.parent=b},i.BoneData.prototype={length:0,x:0,y:0,rotation:0,scaleX:1,scaleY:1},i.SlotData=function(a,b){this.name=a,this.boneData=b},i.SlotData.prototype={r:1,g:1,b:1,a:1,attachmentName:null},i.Bone=function(a,b){this.data=a,this.parent=b,this.setToSetupPose()},i.Bone.yDown=!1,i.Bone.prototype={x:0,y:0,rotation:0,scaleX:1,scaleY:1,m00:0,m01:0,worldX:0,m10:0,m11:0,worldY:0,worldRotation:0,worldScaleX:1,worldScaleY:1,updateWorldTransform:function(a,b){var c=this.parent;null!=c?(this.worldX=this.x*c.m00+this.y*c.m01+c.worldX,this.worldY=this.x*c.m10+this.y*c.m11+c.worldY,this.worldScaleX=c.worldScaleX*this.scaleX,this.worldScaleY=c.worldScaleY*this.scaleY,this.worldRotation=c.worldRotation+this.rotation):(this.worldX=this.x,this.worldY=this.y,this.worldScaleX=this.scaleX,this.worldScaleY=this.scaleY,this.worldRotation=this.rotation);var d=this.worldRotation*Math.PI/180,e=Math.cos(d),f=Math.sin(d);this.m00=e*this.worldScaleX,this.m10=f*this.worldScaleX,this.m01=-f*this.worldScaleY,this.m11=e*this.worldScaleY,a&&(this.m00=-this.m00,this.m01=-this.m01),b&&(this.m10=-this.m10,this.m11=-this.m11),i.Bone.yDown&&(this.m10=-this.m10,this.m11=-this.m11)},setToSetupPose:function(){var a=this.data;this.x=a.x,this.y=a.y,this.rotation=a.rotation,this.scaleX=a.scaleX,this.scaleY=a.scaleY}},i.Slot=function(a,b,c){this.data=a,this.skeleton=b,this.bone=c,this.setToSetupPose()},i.Slot.prototype={r:1,g:1,b:1,a:1,_attachmentTime:0,attachment:null,setAttachment:function(a){this.attachment=a,this._attachmentTime=this.skeleton.time},setAttachmentTime:function(a){this._attachmentTime=this.skeleton.time-a},getAttachmentTime:function(){return this.skeleton.time-this._attachmentTime},setToSetupPose:function(){var a=this.data;this.r=a.r,this.g=a.g,this.b=a.b,this.a=a.a;for(var b=this.skeleton.data.slots,c=0,d=b.length;d>c;c++)if(b[c]==a){this.setAttachment(a.attachmentName?this.skeleton.getAttachmentBySlotIndex(c,a.attachmentName):null);break}}},i.Skin=function(a){this.name=a,this.attachments={}},i.Skin.prototype={addAttachment:function(a,b,c){this.attachments[a+":"+b]=c},getAttachment:function(a,b){return this.attachments[a+":"+b]},_attachAll:function(a,b){for(var c in b.attachments){var d=c.indexOf(":"),e=parseInt(c.substring(0,d),10),f=c.substring(d+1),g=a.slots[e];if(g.attachment&&g.attachment.name==f){var h=this.getAttachment(e,f);h&&g.setAttachment(h)}}}},i.Animation=function(a,b,c){this.name=a,this.timelines=b,this.duration=c},i.Animation.prototype={apply:function(a,b,c){c&&this.duration&&(b%=this.duration);for(var d=this.timelines,e=0,f=d.length;f>e;e++)d[e].apply(a,b,1)},mix:function(a,b,c,d){c&&this.duration&&(b%=this.duration);for(var e=this.timelines,f=0,g=e.length;g>f;f++)e[f].apply(a,b,d)}},i.binarySearch=function(a,b,c){var d=0,e=Math.floor(a.length/c)-2;if(!e)return c;for(var f=e>>>1;;){if(a[(f+1)*c]<=b?d=f+1:e=f,d==e)return(d+1)*c;f=d+e>>>1}},i.linearSearch=function(a,b,c){for(var d=0,e=a.length-c;e>=d;d+=c)if(a[d]>b)return d;return-1},i.Curves=function(a){this.curves=[],this.curves.length=6*(a-1)},i.Curves.prototype={setLinear:function(a){this.curves[6*a]=0},setStepped:function(a){this.curves[6*a]=-1},setCurve:function(a,b,c,d,e){var f=.1,g=f*f,h=g*f,i=3*f,j=3*g,k=6*g,l=6*h,m=2*-b+d,n=2*-c+e,o=3*(b-d)+1,p=3*(c-e)+1,q=6*a,r=this.curves;r[q]=b*i+m*j+o*h,r[q+1]=c*i+n*j+p*h,r[q+2]=m*k+o*l,r[q+3]=n*k+p*l,r[q+4]=o*l,r[q+5]=p*l},getCurvePercent:function(a,b){b=0>b?0:b>1?1:b;var c=6*a,d=this.curves,e=d[c];if(!e)return b;if(-1==e)return 0;for(var f=d[c+1],g=d[c+2],h=d[c+3],i=d[c+4],j=d[c+5],k=e,l=f,m=8;;){if(k>=b){var n=k-e,o=l-f;return o+(l-o)*(b-n)/(k-n)}if(!m)break;m--,e+=g,f+=h,g+=i,h+=j,k+=e,l+=f}return l+(1-l)*(b-k)/(1-k)}},i.RotateTimeline=function(a){this.curves=new i.Curves(a),this.frames=[],this.frames.length=2*a},i.RotateTimeline.prototype={boneIndex:0,getFrameCount:function(){return this.frames.length/2},setFrame:function(a,b,c){a*=2,this.frames[a]=b,this.frames[a+1]=c},apply:function(a,b,c){var d,e=this.frames;if(!(b<e[0])){var f=a.bones[this.boneIndex];if(b>=e[e.length-2]){for(d=f.data.rotation+e[e.length-1]-f.rotation;d>180;)d-=360;for(;-180>d;)d+=360;return f.rotation+=d*c,void 0}var g=i.binarySearch(e,b,2),h=e[g-1],j=e[g],k=1-(b-j)/(e[g-2]-j);for(k=this.curves.getCurvePercent(g/2-1,k),d=e[g+1]-h;d>180;)d-=360;for(;-180>d;)d+=360;for(d=f.data.rotation+(h+d*k)-f.rotation;d>180;)d-=360;for(;-180>d;)d+=360;f.rotation+=d*c}}},i.TranslateTimeline=function(a){this.curves=new i.Curves(a),this.frames=[],this.frames.length=3*a},i.TranslateTimeline.prototype={boneIndex:0,getFrameCount:function(){return this.frames.length/3},setFrame:function(a,b,c,d){a*=3,this.frames[a]=b,this.frames[a+1]=c,this.frames[a+2]=d},apply:function(a,b,c){var d=this.frames;if(!(b<d[0])){var e=a.bones[this.boneIndex];if(b>=d[d.length-3])return e.x+=(e.data.x+d[d.length-2]-e.x)*c,e.y+=(e.data.y+d[d.length-1]-e.y)*c,void 0;var f=i.binarySearch(d,b,3),g=d[f-2],h=d[f-1],j=d[f],k=1-(b-j)/(d[f+-3]-j);k=this.curves.getCurvePercent(f/3-1,k),e.x+=(e.data.x+g+(d[f+1]-g)*k-e.x)*c,e.y+=(e.data.y+h+(d[f+2]-h)*k-e.y)*c}}},i.ScaleTimeline=function(a){this.curves=new i.Curves(a),this.frames=[],this.frames.length=3*a},i.ScaleTimeline.prototype={boneIndex:0,getFrameCount:function(){return this.frames.length/3},setFrame:function(a,b,c,d){a*=3,this.frames[a]=b,this.frames[a+1]=c,this.frames[a+2]=d},apply:function(a,b,c){var d=this.frames;if(!(b<d[0])){var e=a.bones[this.boneIndex];if(b>=d[d.length-3])return e.scaleX+=(e.data.scaleX-1+d[d.length-2]-e.scaleX)*c,e.scaleY+=(e.data.scaleY-1+d[d.length-1]-e.scaleY)*c,void 0;var f=i.binarySearch(d,b,3),g=d[f-2],h=d[f-1],j=d[f],k=1-(b-j)/(d[f+-3]-j);k=this.curves.getCurvePercent(f/3-1,k),e.scaleX+=(e.data.scaleX-1+g+(d[f+1]-g)*k-e.scaleX)*c,e.scaleY+=(e.data.scaleY-1+h+(d[f+2]-h)*k-e.scaleY)*c}}},i.ColorTimeline=function(a){this.curves=new i.Curves(a),this.frames=[],this.frames.length=5*a},i.ColorTimeline.prototype={slotIndex:0,getFrameCount:function(){return this.frames.length/2},setFrame:function(c,d){c*=5,this.frames[c]=d,this.frames[c+1]=r,this.frames[c+2]=g,this.frames[c+3]=b,this.frames[c+4]=a},apply:function(a,b,c){var d=this.frames;if(!(b<d[0])){var e=a.slots[this.slotIndex];if(b>=d[d.length-5]){var f=d.length-1;return e.r=d[f-3],e.g=d[f-2],e.b=d[f-1],e.a=d[f],void 0}var g=i.binarySearch(d,b,5),h=d[g-4],j=d[g-3],k=d[g-2],l=d[g-1],m=d[g],n=1-(b-m)/(d[g-5]-m);n=this.curves.getCurvePercent(g/5-1,n);var o=h+(d[g+1]-h)*n,p=j+(d[g+2]-j)*n,q=k+(d[g+3]-k)*n,r=l+(d[g+4]-l)*n;1>c?(e.r+=(o-e.r)*c,e.g+=(p-e.g)*c,e.b+=(q-e.b)*c,e.a+=(r-e.a)*c):(e.r=o,e.g=p,e.b=q,e.a=r)}}},i.AttachmentTimeline=function(a){this.curves=new i.Curves(a),this.frames=[],this.frames.length=a,this.attachmentNames=[],this.attachmentNames.length=a},i.AttachmentTimeline.prototype={slotIndex:0,getFrameCount:function(){return this.frames.length},setFrame:function(a,b,c){this.frames[a]=b,this.attachmentNames[a]=c},apply:function(a,b){var c=this.frames;if(!(b<c[0])){var d;d=b>=c[c.length-1]?c.length-1:i.binarySearch(c,b,1)-1;var e=this.attachmentNames[d];a.slots[this.slotIndex].setAttachment(e?a.getAttachmentBySlotIndex(this.slotIndex,e):null)}}},i.SkeletonData=function(){this.bones=[],this.slots=[],this.skins=[],this.animations=[]},i.SkeletonData.prototype={defaultSkin:null,findBone:function(a){for(var b=this.bones,c=0,d=b.length;d>c;c++)if(b[c].name==a)return b[c];return null},findBoneIndex:function(a){for(var b=this.bones,c=0,d=b.length;d>c;c++)if(b[c].name==a)return c;return-1},findSlot:function(a){for(var b=this.slots,c=0,d=b.length;d>c;c++)if(b[c].name==a)return slot[c];return null},findSlotIndex:function(a){for(var b=this.slots,c=0,d=b.length;d>c;c++)if(b[c].name==a)return c;return-1},findSkin:function(a){for(var b=this.skins,c=0,d=b.length;d>c;c++)if(b[c].name==a)return b[c];return null},findAnimation:function(a){for(var b=this.animations,c=0,d=b.length;d>c;c++)if(b[c].name==a)return b[c];return null}},i.Skeleton=function(a){this.data=a,this.bones=[];for(var b=0,c=a.bones.length;c>b;b++){var d=a.bones[b],e=d.parent?this.bones[a.bones.indexOf(d.parent)]:null;this.bones.push(new i.Bone(d,e))}for(this.slots=[],this.drawOrder=[],b=0,c=a.slots.length;c>b;b++){var f=a.slots[b],g=this.bones[a.bones.indexOf(f.boneData)],h=new i.Slot(f,this,g);this.slots.push(h),this.drawOrder.push(h)}},i.Skeleton.prototype={x:0,y:0,skin:null,r:1,g:1,b:1,a:1,time:0,flipX:!1,flipY:!1,updateWorldTransform:function(){for(var a=this.flipX,b=this.flipY,c=this.bones,d=0,e=c.length;e>d;d++)c[d].updateWorldTransform(a,b)},setToSetupPose:function(){this.setBonesToSetupPose(),this.setSlotsToSetupPose()},setBonesToSetupPose:function(){for(var a=this.bones,b=0,c=a.length;c>b;b++)a[b].setToSetupPose()},setSlotsToSetupPose:function(){for(var a=this.slots,b=0,c=a.length;c>b;b++)a[b].setToSetupPose(b)},getRootBone:function(){return this.bones.length?this.bones[0]:null},findBone:function(a){for(var b=this.bones,c=0,d=b.length;d>c;c++)if(b[c].data.name==a)return b[c];return null},findBoneIndex:function(a){for(var b=this.bones,c=0,d=b.length;d>c;c++)if(b[c].data.name==a)return c;return-1},findSlot:function(a){for(var b=this.slots,c=0,d=b.length;d>c;c++)if(b[c].data.name==a)return b[c];return null},findSlotIndex:function(a){for(var b=this.slots,c=0,d=b.length;d>c;c++)if(b[c].data.name==a)return c;return-1},setSkinByName:function(a){var b=this.data.findSkin(a);if(!b)throw"Skin not found: "+a;this.setSkin(b)},setSkin:function(a){this.skin&&a&&a._attachAll(this,this.skin),this.skin=a},getAttachmentBySlotName:function(a,b){return this.getAttachmentBySlotIndex(this.data.findSlotIndex(a),b)},getAttachmentBySlotIndex:function(a,b){if(this.skin){var c=this.skin.getAttachment(a,b);if(c)return c}return this.data.defaultSkin?this.data.defaultSkin.getAttachment(a,b):null},setAttachment:function(a,b){for(var c=this.slots,d=0,e=c.size;e>d;d++){var f=c[d];if(f.data.name==a){var g=null;if(b&&(g=this.getAttachment(d,b),null==g))throw"Attachment not found: "+b+", for slot: "+a;return f.setAttachment(g),void 0}}throw"Slot not found: "+a},update:function(a){time+=a}},i.AttachmentType={region:0},i.RegionAttachment=function(){this.offset=[],this.offset.length=8,this.uvs=[],this.uvs.length=8},i.RegionAttachment.prototype={x:0,y:0,rotation:0,scaleX:1,scaleY:1,width:0,height:0,rendererObject:null,regionOffsetX:0,regionOffsetY:0,regionWidth:0,regionHeight:0,regionOriginalWidth:0,regionOriginalHeight:0,setUVs:function(a,b,c,d,e){var f=this.uvs;e?(f[2]=a,f[3]=d,f[4]=a,f[5]=b,f[6]=c,f[7]=b,f[0]=c,f[1]=d):(f[0]=a,f[1]=d,f[2]=a,f[3]=b,f[4]=c,f[5]=b,f[6]=c,f[7]=d)},updateOffset:function(){var a=this.width/this.regionOriginalWidth*this.scaleX,b=this.height/this.regionOriginalHeight*this.scaleY,c=-this.width/2*this.scaleX+this.regionOffsetX*a,d=-this.height/2*this.scaleY+this.regionOffsetY*b,e=c+this.regionWidth*a,f=d+this.regionHeight*b,g=this.rotation*Math.PI/180,h=Math.cos(g),i=Math.sin(g),j=c*h+this.x,k=c*i,l=d*h+this.y,m=d*i,n=e*h+this.x,o=e*i,p=f*h+this.y,q=f*i,r=this.offset;r[0]=j-m,r[1]=l+k,r[2]=j-q,r[3]=p+k,r[4]=n-q,r[5]=p+o,r[6]=n-m,r[7]=l+o},computeVertices:function(a,b,c,d){a+=c.worldX,b+=c.worldY;var e=c.m00,f=c.m01,g=c.m10,h=c.m11,i=this.offset;d[0]=i[0]*e+i[1]*f+a,d[1]=i[0]*g+i[1]*h+b,d[2]=i[2]*e+i[3]*f+a,d[3]=i[2]*g+i[3]*h+b,d[4]=i[4]*e+i[5]*f+a,d[5]=i[4]*g+i[5]*h+b,d[6]=i[6]*e+i[7]*f+a,d[7]=i[6]*g+i[7]*h+b}},i.AnimationStateData=function(a){this.skeletonData=a,this.animationToMixTime={}},i.AnimationStateData.prototype={defaultMix:0,setMixByName:function(a,b,c){var d=this.skeletonData.findAnimation(a);if(!d)throw"Animation not found: "+a;var e=this.skeletonData.findAnimation(b);if(!e)throw"Animation not found: "+b;this.setMix(d,e,c)},setMix:function(a,b,c){this.animationToMixTime[a.name+":"+b.name]=c},getMix:function(a,b){var c=this.animationToMixTime[a.name+":"+b.name];return c?c:this.defaultMix}},i.AnimationState=function(a){this.data=a,this.queue=[]},i.AnimationState.prototype={current:null,previous:null,currentTime:0,previousTime:0,currentLoop:!1,previousLoop:!1,mixTime:0,mixDuration:0,update:function(a){if(this.currentTime+=a,this.previousTime+=a,this.mixTime+=a,this.queue.length>0){var b=this.queue[0];this.currentTime>=b.delay&&(this._setAnimation(b.animation,b.loop),this.queue.shift())}},apply:function(a){if(this.current)if(this.previous){this.previous.apply(a,this.previousTime,this.previousLoop);var b=this.mixTime/this.mixDuration;b>=1&&(b=1,this.previous=null),this.current.mix(a,this.currentTime,this.currentLoop,b)}else this.current.apply(a,this.currentTime,this.currentLoop)},clearAnimation:function(){this.previous=null,this.current=null,this.queue.length=0},_setAnimation:function(a,b){this.previous=null,a&&this.current&&(this.mixDuration=this.data.getMix(this.current,a),this.mixDuration>0&&(this.mixTime=0,this.previous=this.current,this.previousTime=this.currentTime,this.previousLoop=this.currentLoop)),this.current=a,this.currentLoop=b,this.currentTime=0},setAnimationByName:function(a,b){var c=this.data.skeletonData.findAnimation(a);if(!c)throw"Animation not found: "+a;this.setAnimation(c,b)},setAnimation:function(a,b){this.queue.length=0,this._setAnimation(a,b)},addAnimationByName:function(a,b,c){var d=this.data.skeletonData.findAnimation(a);if(!d)throw"Animation not found: "+a;this.addAnimation(d,b,c)},addAnimation:function(a,b,c){var d={};if(d.animation=a,d.loop=b,!c||0>=c){var e=this.queue.length?this.queue[this.queue.length-1].animation:this.current;c=null!=e?e.duration-this.data.getMix(e,a)+(c||0):0}d.delay=c,this.queue.push(d)},isComplete:function(){return!this.current||this.currentTime>=this.current.duration}},i.SkeletonJson=function(a){this.attachmentLoader=a},i.SkeletonJson.prototype={scale:1,readSkeletonData:function(a){for(var b,c=new i.SkeletonData,d=a.bones,e=0,f=d.length;f>e;e++){var g=d[e],h=null;if(g.parent&&(h=c.findBone(g.parent),!h))throw"Parent bone not found: "+g.parent;b=new i.BoneData(g.name,h),b.length=(g.length||0)*this.scale,b.x=(g.x||0)*this.scale,b.y=(g.y||0)*this.scale,b.rotation=g.rotation||0,b.scaleX=g.scaleX||1,b.scaleY=g.scaleY||1,c.bones.push(b)}var j=a.slots;for(e=0,f=j.length;f>e;e++){var k=j[e];if(b=c.findBone(k.bone),!b)throw"Slot bone not found: "+k.bone;var l=new i.SlotData(k.name,b),m=k.color;m&&(l.r=i.SkeletonJson.toColor(m,0),l.g=i.SkeletonJson.toColor(m,1),l.b=i.SkeletonJson.toColor(m,2),l.a=i.SkeletonJson.toColor(m,3)),l.attachmentName=k.attachment,c.slots.push(l)}var n=a.skins;for(var o in n)if(n.hasOwnProperty(o)){var p=n[o],q=new i.Skin(o);for(var r in p)if(p.hasOwnProperty(r)){var s=c.findSlotIndex(r),t=p[r];for(var u in t)if(t.hasOwnProperty(u)){var v=this.readAttachment(q,u,t[u]);null!=v&&q.addAttachment(s,u,v)}}c.skins.push(q),"default"==q.name&&(c.defaultSkin=q)}var w=a.animations;for(var x in w)w.hasOwnProperty(x)&&this.readAnimation(x,w[x],c);return c},readAttachment:function(a,b,c){b=c.name||b;var d=i.AttachmentType[c.type||"region"];if(d==i.AttachmentType.region){var e=new i.RegionAttachment;return e.x=(c.x||0)*this.scale,e.y=(c.y||0)*this.scale,e.scaleX=c.scaleX||1,e.scaleY=c.scaleY||1,e.rotation=c.rotation||0,e.width=(c.width||32)*this.scale,e.height=(c.height||32)*this.scale,e.updateOffset(),e.rendererObject={},e.rendererObject.name=b,e.rendererObject.scale={},e.rendererObject.scale.x=e.scaleX,e.rendererObject.scale.y=e.scaleY,e.rendererObject.rotation=-e.rotation*Math.PI/180,e}throw"Unknown attachment type: "+d},readAnimation:function(a,b,c){var d,e,f,g,h,j,k,l=[],m=0,n=b.bones;for(var o in n)if(n.hasOwnProperty(o)){var p=c.findBoneIndex(o);if(-1==p)throw"Bone not found: "+o;var q=n[o];for(f in q)if(q.hasOwnProperty(f))if(h=q[f],"rotate"==f){for(e=new i.RotateTimeline(h.length),e.boneIndex=p,d=0,j=0,k=h.length;k>j;j++)g=h[j],e.setFrame(d,g.time,g.angle),i.SkeletonJson.readCurve(e,d,g),d++;l.push(e),m=Math.max(m,e.frames[2*e.getFrameCount()-2])}else{if("translate"!=f&&"scale"!=f)throw"Invalid timeline type for a bone: "+f+" ("+o+")";var r=1;for("scale"==f?e=new i.ScaleTimeline(h.length):(e=new i.TranslateTimeline(h.length),r=this.scale),e.boneIndex=p,d=0,j=0,k=h.length;k>j;j++){g=h[j];var s=(g.x||0)*r,t=(g.y||0)*r;e.setFrame(d,g.time,s,t),i.SkeletonJson.readCurve(e,d,g),d++}l.push(e),m=Math.max(m,e.frames[3*e.getFrameCount()-3])}}var u=b.slots;for(var v in u)if(u.hasOwnProperty(v)){var w=u[v],x=c.findSlotIndex(v);for(f in w)if(w.hasOwnProperty(f))if(h=w[f],"color"==f){for(e=new i.ColorTimeline(h.length),e.slotIndex=x,d=0,j=0,k=h.length;k>j;j++){g=h[j];var y=g.color,z=i.SkeletonJson.toColor(y,0),A=i.SkeletonJson.toColor(y,1),B=i.SkeletonJson.toColor(y,2),C=i.SkeletonJson.toColor(y,3);e.setFrame(d,g.time,z,A,B,C),i.SkeletonJson.readCurve(e,d,g),d++}l.push(e),m=Math.max(m,e.frames[5*e.getFrameCount()-5])}else{if("attachment"!=f)throw"Invalid timeline type for a slot: "+f+" ("+v+")";for(e=new i.AttachmentTimeline(h.length),e.slotIndex=x,d=0,j=0,k=h.length;k>j;j++)g=h[j],e.setFrame(d++,g.time,g.name);l.push(e),m=Math.max(m,e.frames[e.getFrameCount()-1])}}c.animations.push(new i.Animation(a,l,m))}},i.SkeletonJson.readCurve=function(a,b,c){var d=c.curve;d&&("stepped"==d?a.curves.setStepped(b):d instanceof Array&&a.curves.setCurve(b,d[0],d[1],d[2],d[3]))},i.SkeletonJson.toColor=function(a,b){if(8!=a.length)throw"Color hexidecimal length must be 8, recieved: "+a;return parseInt(a.substring(2*b,2),16)/255},i.Atlas=function(a,b){this.textureLoader=b,this.pages=[],this.regions=[];var c=new i.AtlasReader(a),d=[];d.length=4;for(var e=null;;){var f=c.readLine();if(null==f)break;if(f=c.trim(f),f.length)if(e){var g=new i.AtlasRegion;g.name=f,g.page=e,g.rotate="true"==c.readValue(),c.readTuple(d);var h=parseInt(d[0],10),j=parseInt(d[1],10);c.readTuple(d);var k=parseInt(d[0],10),l=parseInt(d[1],10);g.u=h/e.width,g.v=j/e.height,g.rotate?(g.u2=(h+l)/e.width,g.v2=(j+k)/e.height):(g.u2=(h+k)/e.width,g.v2=(j+l)/e.height),g.x=h,g.y=j,g.width=Math.abs(k),g.height=Math.abs(l),4==c.readTuple(d)&&(g.splits=[parseInt(d[0],10),parseInt(d[1],10),parseInt(d[2],10),parseInt(d[3],10)],4==c.readTuple(d)&&(g.pads=[parseInt(d[0],10),parseInt(d[1],10),parseInt(d[2],10),parseInt(d[3],10)],c.readTuple(d))),g.originalWidth=parseInt(d[0],10),g.originalHeight=parseInt(d[1],10),c.readTuple(d),g.offsetX=parseInt(d[0],10),g.offsetY=parseInt(d[1],10),g.index=parseInt(c.readValue(),10),this.regions.push(g)}else{e=new i.AtlasPage,e.name=f,e.format=i.Atlas.Format[c.readValue()],c.readTuple(d),e.minFilter=i.Atlas.TextureFilter[d[0]],e.magFilter=i.Atlas.TextureFilter[d[1]];var m=c.readValue();e.uWrap=i.Atlas.TextureWrap.clampToEdge,e.vWrap=i.Atlas.TextureWrap.clampToEdge,"x"==m?e.uWrap=i.Atlas.TextureWrap.repeat:"y"==m?e.vWrap=i.Atlas.TextureWrap.repeat:"xy"==m&&(e.uWrap=e.vWrap=i.Atlas.TextureWrap.repeat),b.load(e,f),this.pages.push(e)}else e=null}},i.Atlas.prototype={findRegion:function(a){for(var b=this.regions,c=0,d=b.length;d>c;c++)if(b[c].name==a)return b[c];return null},dispose:function(){for(var a=this.pages,b=0,c=a.length;c>b;b++)this.textureLoader.unload(a[b].rendererObject)},updateUVs:function(a){for(var b=this.regions,c=0,d=b.length;d>c;c++){var e=b[c];e.page==a&&(e.u=e.x/a.width,e.v=e.y/a.height,e.rotate?(e.u2=(e.x+e.height)/a.width,e.v2=(e.y+e.width)/a.height):(e.u2=(e.x+e.width)/a.width,e.v2=(e.y+e.height)/a.height))}}},i.Atlas.Format={alpha:0,intensity:1,luminanceAlpha:2,rgb565:3,rgba4444:4,rgb888:5,rgba8888:6},i.Atlas.TextureFilter={nearest:0,linear:1,mipMap:2,mipMapNearestNearest:3,mipMapLinearNearest:4,mipMapNearestLinear:5,mipMapLinearLinear:6},i.Atlas.TextureWrap={mirroredRepeat:0,clampToEdge:1,repeat:2},i.AtlasPage=function(){},i.AtlasPage.prototype={name:null,format:null,minFilter:null,magFilter:null,uWrap:null,vWrap:null,rendererObject:null,width:0,height:0},i.AtlasRegion=function(){},i.AtlasRegion.prototype={page:null,name:null,x:0,y:0,width:0,height:0,u:0,v:0,u2:0,v2:0,offsetX:0,offsetY:0,originalWidth:0,originalHeight:0,index:0,rotate:!1,splits:null,pads:null},i.AtlasReader=function(a){this.lines=a.split(/\r\n|\r|\n/)},i.AtlasReader.prototype={index:0,trim:function(a){return a.replace(/^\s+|\s+$/g,"")},readLine:function(){return this.index>=this.lines.length?null:this.lines[this.index++]},readValue:function(){var a=this.readLine(),b=a.indexOf(":");if(-1==b)throw"Invalid line: "+a;return this.trim(a.substring(b+1))},readTuple:function(a){var b=this.readLine(),c=b.indexOf(":");if(-1==c)throw"Invalid line: "+b;for(var d=0,e=c+1;3>d;d++){var f=b.indexOf(",",e);if(-1==f){if(!d)throw"Invalid line: "+b;break}a[d]=this.trim(b.substr(e,f-e)),e=f+1}return a[d]=this.trim(b.substring(e)),d+1}},i.AtlasAttachmentLoader=function(a){this.atlas=a},i.AtlasAttachmentLoader.prototype={newAttachment:function(a,b,c){switch(b){case i.AttachmentType.region:var d=this.atlas.findRegion(c);if(!d)throw"Region not found in atlas: "+c+" ("+b+")";var e=new i.RegionAttachment(c);return e.rendererObject=d,e.setUVs(d.u,d.v,d.u2,d.v2,d.rotate),e.regionOffsetX=d.offsetX,e.regionOffsetY=d.offsetY,e.regionWidth=d.width,e.regionHeight=d.height,e.regionOriginalWidth=d.originalWidth,e.regionOriginalHeight=d.originalHeight,e}throw"Unknown attachment type: "+b}},i.Bone.yDown=!0,d.AnimCache={},d.Spine=function(a){if(d.DisplayObjectContainer.call(this),this.spineData=d.AnimCache[a],!this.spineData)throw new Error("Spine data must be preloaded using PIXI.SpineLoader or PIXI.AssetLoader: "+a);this.skeleton=new i.Skeleton(this.spineData),this.skeleton.updateWorldTransform(),this.stateData=new i.AnimationStateData(this.spineData),this.state=new i.AnimationState(this.stateData),this.slotContainers=[];for(var b=0,c=this.skeleton.drawOrder.length;c>b;b++){var e=this.skeleton.drawOrder[b],f=e.attachment,g=new d.DisplayObjectContainer;if(this.slotContainers.push(g),this.addChild(g),f instanceof i.RegionAttachment){var h=f.rendererObject.name,j=this.createSprite(e,f.rendererObject);e.currentSprite=j,e.currentSpriteName=h,g.addChild(j)}}},d.Spine.prototype=Object.create(d.DisplayObjectContainer.prototype),d.Spine.prototype.constructor=d.Spine,d.Spine.prototype.updateTransform=function(){this.lastTime=this.lastTime||Date.now();var a=.001*(Date.now()-this.lastTime);this.lastTime=Date.now(),this.state.update(a),this.state.apply(this.skeleton),this.skeleton.updateWorldTransform();for(var b=this.skeleton.drawOrder,c=0,e=b.length;e>c;c++){var f=b[c],g=f.attachment,h=this.slotContainers[c];if(g instanceof i.RegionAttachment){if(g.rendererObject&&(!f.currentSpriteName||f.currentSpriteName!=g.name)){var j=g.rendererObject.name;if(void 0!==f.currentSprite&&(f.currentSprite.visible=!1),f.sprites=f.sprites||{},void 0!==f.sprites[j])f.sprites[j].visible=!0;else{var k=this.createSprite(f,g.rendererObject);h.addChild(k)}f.currentSprite=f.sprites[j],f.currentSpriteName=j}h.visible=!0;var l=f.bone;h.position.x=l.worldX+g.x*l.m00+g.y*l.m01,h.position.y=l.worldY+g.x*l.m10+g.y*l.m11,h.scale.x=l.worldScaleX,h.scale.y=l.worldScaleY,h.rotation=-(f.bone.worldRotation*Math.PI/180)}else h.visible=!1}d.DisplayObjectContainer.prototype.updateTransform.call(this)},d.Spine.prototype.createSprite=function(a,b){var c=d.TextureCache[b.name]?b.name:b.name+".png",e=new d.Sprite(d.Texture.fromFrame(c));return e.scale=b.scale,e.rotation=b.rotation,e.anchor.x=e.anchor.y=.5,a.sprites=a.sprites||{},a.sprites[b.name]=e,e},d.BaseTextureCache={},d.texturesToUpdate=[],d.texturesToDestroy=[],d.BaseTextureCacheIdGenerator=0,d.BaseTexture=function(a,b){if(d.EventTarget.call(this),this.width=100,this.height=100,this.scaleMode=b||d.scaleModes.DEFAULT,this.hasLoaded=!1,this.source=a,a){if(this.source.complete||this.source.getContext)this.hasLoaded=!0,this.width=this.source.width,this.height=this.source.height,d.texturesToUpdate.push(this);else{var c=this;this.source.onload=function(){c.hasLoaded=!0,c.width=c.source.width,c.height=c.source.height,d.texturesToUpdate.push(c),c.dispatchEvent({type:"loaded",content:c})}}this.imageUrl=null,this._powerOf2=!1,this.id=d.BaseTextureCacheIdGenerator++,this._glTextures=[]}},d.BaseTexture.prototype.constructor=d.BaseTexture,d.BaseTexture.prototype.destroy=function(){this.imageUrl&&(delete d.BaseTextureCache[this.imageUrl],this.imageUrl=null,this.source.src=null),this.source=null,d.texturesToDestroy.push(this)},d.BaseTexture.prototype.updateSourceImage=function(a){this.hasLoaded=!1,this.source.src=null,this.source.src=a},d.BaseTexture.fromImage=function(a,b,c){var e=d.BaseTextureCache[a];if(b=!b,!e){var f=new Image;b&&(f.crossOrigin=""),f.src=a,e=new d.BaseTexture(f,c),e.imageUrl=a,d.BaseTextureCache[a]=e}return e},d.BaseTexture.fromCanvas=function(a,b){a._pixiId||(a._pixiId="canvas_"+d.TextureCacheIdGenerator++);var c=d.BaseTextureCache[a._pixiId];return c||(c=new d.BaseTexture(a,b),d.BaseTextureCache[a._pixiId]=c),c},d.TextureCache={},d.FrameCache={},d.TextureCacheIdGenerator=0,d.Texture=function(a,b){if(d.EventTarget.call(this),b||(this.noFrame=!0,b=new d.Rectangle(0,0,1,1)),a instanceof d.Texture&&(a=a.baseTexture),this.baseTexture=a,this.frame=b,this.trim=null,this.scope=this,a.hasLoaded)this.noFrame&&(b=new d.Rectangle(0,0,a.width,a.height)),this.setFrame(b);else{var c=this;a.addEventListener("loaded",function(){c.onBaseTextureLoaded()})}},d.Texture.prototype.constructor=d.Texture,d.Texture.prototype.onBaseTextureLoaded=function(){var a=this.baseTexture;a.removeEventListener("loaded",this.onLoaded),this.noFrame&&(this.frame=new d.Rectangle(0,0,a.width,a.height)),this.setFrame(this.frame),this.scope.dispatchEvent({type:"update",content:this})},d.Texture.prototype.destroy=function(a){a&&this.baseTexture.destroy()},d.Texture.prototype.setFrame=function(a){if(this.frame=a,this.width=a.width,this.height=a.height,a.x+a.width>this.baseTexture.width||a.y+a.height>this.baseTexture.height)throw new Error("Texture Error: frame does not fit inside the base Texture dimensions "+this);this.updateFrame=!0,d.Texture.frameUpdates.push(this)},d.Texture.prototype._updateWebGLuvs=function(){this._uvs||(this._uvs=new d.TextureUvs);var a=this.frame,b=this.baseTexture.width,c=this.baseTexture.height;this._uvs.x0=a.x/b,this._uvs.y0=a.y/c,this._uvs.x1=(a.x+a.width)/b,this._uvs.y1=a.y/c,this._uvs.x2=(a.x+a.width)/b,this._uvs.y2=(a.y+a.height)/c,this._uvs.x3=a.x/b,this._uvs.y3=(a.y+a.height)/c},d.Texture.fromImage=function(a,b,c){var e=d.TextureCache[a];return e||(e=new d.Texture(d.BaseTexture.fromImage(a,b,c)),d.TextureCache[a]=e),e},d.Texture.fromFrame=function(a){var b=d.TextureCache[a];if(!b)throw new Error('The frameId "'+a+'" does not exist in the texture cache ');return b},d.Texture.fromCanvas=function(a,b){var c=d.BaseTexture.fromCanvas(a,b);return new d.Texture(c)},d.Texture.addTextureToCache=function(a,b){d.TextureCache[b]=a},d.Texture.removeTextureFromCache=function(a){var b=d.TextureCache[a];return d.TextureCache[a]=null,b},d.Texture.frameUpdates=[],d.TextureUvs=function(){this.x0=0,this.y0=0,this.x1=0,this.y1=0,this.x2=0,this.y2=0,this.x3=0,this.y4=0},d.RenderTexture=function(a,b,c){if(d.EventTarget.call(this),this.width=a||100,this.height=b||100,this.frame=new d.Rectangle(0,0,this.width,this.height),this.baseTexture=new d.BaseTexture,this.baseTexture.width=this.width,this.baseTexture.height=this.height,this.baseTexture._glTextures=[],this.baseTexture.hasLoaded=!0,this.renderer=c||d.defaultRenderer,this.renderer.type===d.WEBGL_RENDERER){var e=this.renderer.gl;this.textureBuffer=new d.FilterTexture(e,this.width,this.height),this.baseTexture._glTextures[e.id]=this.textureBuffer.texture,this.render=this.renderWebGL,this.projection=new d.Point(this.width/2,-this.height/2)}else this.render=this.renderCanvas,this.textureBuffer=new d.CanvasBuffer(this.width,this.height),this.baseTexture.source=this.textureBuffer.canvas;d.Texture.frameUpdates.push(this)},d.RenderTexture.prototype=Object.create(d.Texture.prototype),d.RenderTexture.prototype.constructor=d.RenderTexture,d.RenderTexture.prototype.resize=function(a,b){if(this.width=a,this.height=b,this.frame.width=this.width,this.frame.height=this.height,this.renderer.type===d.WEBGL_RENDERER){this.projection.x=this.width/2,this.projection.y=-this.height/2;var c=this.renderer.gl;c.bindTexture(c.TEXTURE_2D,this.baseTexture._glTextures[c.id]),c.texImage2D(c.TEXTURE_2D,0,c.RGBA,this.width,this.height,0,c.RGBA,c.UNSIGNED_BYTE,null)}else this.textureBuffer.resize(this.width,this.height);d.Texture.frameUpdates.push(this)},d.RenderTexture.prototype.renderWebGL=function(a,b,c){var e=this.renderer.gl;e.colorMask(!0,!0,!0,!0),e.viewport(0,0,this.width,this.height),e.bindFramebuffer(e.FRAMEBUFFER,this.textureBuffer.frameBuffer),c&&this.textureBuffer.clear();var f=a.children,g=a.worldTransform;a.worldTransform=d.RenderTexture.tempMatrix,a.worldTransform.d=-1,a.worldTransform.ty=-2*this.projection.y,b&&(a.worldTransform.tx=b.x,a.worldTransform.ty-=b.y);for(var h=0,i=f.length;i>h;h++)f[h].updateTransform();d.WebGLRenderer.updateTextures(),this.renderer.renderDisplayObject(a,this.projection,this.textureBuffer.frameBuffer),a.worldTransform=g},d.RenderTexture.prototype.renderCanvas=function(a,b,c){var e=a.children;
a.worldTransform=d.RenderTexture.tempMatrix,b&&(a.worldTransform.tx=b.x,a.worldTransform.ty=b.y);for(var f=0,g=e.length;g>f;f++)e[f].updateTransform();c&&this.textureBuffer.clear();var h=this.textureBuffer.context;this.renderer.renderDisplayObject(a,h),h.setTransform(1,0,0,1,0,0)},d.RenderTexture.tempMatrix=new d.Matrix,d.AssetLoader=function(a,b){d.EventTarget.call(this),this.assetURLs=a,this.crossorigin=b,this.loadersByType={jpg:d.ImageLoader,jpeg:d.ImageLoader,png:d.ImageLoader,gif:d.ImageLoader,json:d.JsonLoader,atlas:d.AtlasLoader,anim:d.SpineLoader,xml:d.BitmapFontLoader,fnt:d.BitmapFontLoader}},d.AssetLoader.prototype.constructor=d.AssetLoader,d.AssetLoader.prototype._getDataType=function(a){var b="data:",c=a.slice(0,b.length).toLowerCase();if(c===b){var d=a.slice(b.length),e=d.indexOf(",");if(-1===e)return null;var f=d.slice(0,e).split(";")[0];return f&&"text/plain"!==f.toLowerCase()?f.split("/").pop().toLowerCase():"txt"}return null},d.AssetLoader.prototype.load=function(){function a(a){b.onAssetLoaded(a.loader)}var b=this;this.loadCount=this.assetURLs.length;for(var c=0;c<this.assetURLs.length;c++){var d=this.assetURLs[c],e=this._getDataType(d);e||(e=d.split("?").shift().split(".").pop().toLowerCase());var f=this.loadersByType[e];if(!f)throw new Error(e+" is an unsupported file type");var g=new f(d,this.crossorigin);g.addEventListener("loaded",a),g.load()}},d.AssetLoader.prototype.onAssetLoaded=function(a){this.loadCount--,this.dispatchEvent({type:"onProgress",content:this,loader:a}),this.onProgress&&this.onProgress(a),this.loadCount||(this.dispatchEvent({type:"onComplete",content:this}),this.onComplete&&this.onComplete())},d.JsonLoader=function(a,b){d.EventTarget.call(this),this.url=a,this.crossorigin=b,this.baseUrl=a.replace(/[^\/]*$/,""),this.loaded=!1},d.JsonLoader.prototype.constructor=d.JsonLoader,d.JsonLoader.prototype.load=function(){this.ajaxRequest=new d.AjaxRequest(this.crossorigin);var a=this;this.ajaxRequest.onreadystatechange=function(){a.onJSONLoaded()},this.ajaxRequest.open("GET",this.url,!0),this.ajaxRequest.overrideMimeType&&this.ajaxRequest.overrideMimeType("application/json"),this.ajaxRequest.send(null)},d.JsonLoader.prototype.onJSONLoaded=function(){if(4===this.ajaxRequest.readyState)if(200===this.ajaxRequest.status||-1===window.location.protocol.indexOf("http"))if(this.json=JSON.parse(this.ajaxRequest.responseText),this.json.frames){var a=this,b=this.baseUrl+this.json.meta.image,c=new d.ImageLoader(b,this.crossorigin),e=this.json.frames;this.texture=c.texture.baseTexture,c.addEventListener("loaded",function(){a.onLoaded()});for(var f in e){var g=e[f].frame;if(g&&(d.TextureCache[f]=new d.Texture(this.texture,{x:g.x,y:g.y,width:g.w,height:g.h}),e[f].trimmed)){var h=d.TextureCache[f],j=e[f].sourceSize,k=e[f].spriteSourceSize;h.trim=new d.Rectangle(k.x,k.y,j.w,j.h)}}c.load()}else if(this.json.bones){var l=new i.SkeletonJson,m=l.readSkeletonData(this.json);d.AnimCache[this.url]=m,this.onLoaded()}else this.onLoaded();else this.onError()},d.JsonLoader.prototype.onLoaded=function(){this.loaded=!0,this.dispatchEvent({type:"loaded",content:this})},d.JsonLoader.prototype.onError=function(){this.dispatchEvent({type:"error",content:this})},d.AtlasLoader=function(a,b){d.EventTarget.call(this),this.url=a,this.baseUrl=a.replace(/[^\/]*$/,""),this.crossorigin=b,this.loaded=!1},d.AtlasLoader.constructor=d.AtlasLoader,d.AtlasLoader.prototype.load=function(){this.ajaxRequest=new d.AjaxRequest,this.ajaxRequest.onreadystatechange=this.onAtlasLoaded.bind(this),this.ajaxRequest.open("GET",this.url,!0),this.ajaxRequest.overrideMimeType&&this.ajaxRequest.overrideMimeType("application/json"),this.ajaxRequest.send(null)},d.AtlasLoader.prototype.onAtlasLoaded=function(){if(4===this.ajaxRequest.readyState)if(200===this.ajaxRequest.status||-1===window.location.href.indexOf("http")){this.atlas={meta:{image:[]},frames:[]};var a=this.ajaxRequest.responseText.split(/\r?\n/),b=-3,c=0,e=null,f=!1,g=0,h=0,i=this.onLoaded.bind(this);for(g=0;g<a.length;g++)if(a[g]=a[g].replace(/^\s+|\s+$/g,""),""===a[g]&&(f=g+1),a[g].length>0){if(f===g)this.atlas.meta.image.push(a[g]),c=this.atlas.meta.image.length-1,this.atlas.frames.push({}),b=-3;else if(b>0)if(b%7===1)null!=e&&(this.atlas.frames[c][e.name]=e),e={name:a[g],frame:{}};else{var j=a[g].split(" ");if(b%7===3)e.frame.x=Number(j[1].replace(",","")),e.frame.y=Number(j[2]);else if(b%7===4)e.frame.w=Number(j[1].replace(",","")),e.frame.h=Number(j[2]);else if(b%7===5){var k={x:0,y:0,w:Number(j[1].replace(",","")),h:Number(j[2])};k.w>e.frame.w||k.h>e.frame.h?(e.trimmed=!0,e.realSize=k):e.trimmed=!1}}b++}if(null!=e&&(this.atlas.frames[c][e.name]=e),this.atlas.meta.image.length>0){for(this.images=[],h=0;h<this.atlas.meta.image.length;h++){var l=this.baseUrl+this.atlas.meta.image[h],m=this.atlas.frames[h];this.images.push(new d.ImageLoader(l,this.crossorigin));for(g in m){var n=m[g].frame;n&&(d.TextureCache[g]=new d.Texture(this.images[h].texture.baseTexture,{x:n.x,y:n.y,width:n.w,height:n.h}),m[g].trimmed&&(d.TextureCache[g].realSize=m[g].realSize,d.TextureCache[g].trim.x=0,d.TextureCache[g].trim.y=0))}}for(this.currentImageId=0,h=0;h<this.images.length;h++)this.images[h].addEventListener("loaded",i);this.images[this.currentImageId].load()}else this.onLoaded()}else this.onError()},d.AtlasLoader.prototype.onLoaded=function(){this.images.length-1>this.currentImageId?(this.currentImageId++,this.images[this.currentImageId].load()):(this.loaded=!0,this.dispatchEvent({type:"loaded",content:this}))},d.AtlasLoader.prototype.onError=function(){this.dispatchEvent({type:"error",content:this})},d.SpriteSheetLoader=function(a,b){d.EventTarget.call(this),this.url=a,this.crossorigin=b,this.baseUrl=a.replace(/[^\/]*$/,""),this.texture=null,this.frames={}},d.SpriteSheetLoader.prototype.constructor=d.SpriteSheetLoader,d.SpriteSheetLoader.prototype.load=function(){var a=this,b=new d.JsonLoader(this.url,this.crossorigin);b.addEventListener("loaded",function(b){a.json=b.content.json,a.onLoaded()}),b.load()},d.SpriteSheetLoader.prototype.onLoaded=function(){this.dispatchEvent({type:"loaded",content:this})},d.ImageLoader=function(a,b){d.EventTarget.call(this),this.texture=d.Texture.fromImage(a,b),this.frames=[]},d.ImageLoader.prototype.constructor=d.ImageLoader,d.ImageLoader.prototype.load=function(){if(this.texture.baseTexture.hasLoaded)this.onLoaded();else{var a=this;this.texture.baseTexture.addEventListener("loaded",function(){a.onLoaded()})}},d.ImageLoader.prototype.onLoaded=function(){this.dispatchEvent({type:"loaded",content:this})},d.ImageLoader.prototype.loadFramedSpriteSheet=function(a,b,c){this.frames=[];for(var e=Math.floor(this.texture.width/a),f=Math.floor(this.texture.height/b),g=0,h=0;f>h;h++)for(var i=0;e>i;i++,g++){var j=new d.Texture(this.texture,{x:i*a,y:h*b,width:a,height:b});this.frames.push(j),c&&(d.TextureCache[c+"-"+g]=j)}if(this.texture.baseTexture.hasLoaded)this.onLoaded();else{var k=this;this.texture.baseTexture.addEventListener("loaded",function(){k.onLoaded()})}},d.BitmapFontLoader=function(a,b){d.EventTarget.call(this),this.url=a,this.crossorigin=b,this.baseUrl=a.replace(/[^\/]*$/,""),this.texture=null},d.BitmapFontLoader.prototype.constructor=d.BitmapFontLoader,d.BitmapFontLoader.prototype.load=function(){this.ajaxRequest=new d.AjaxRequest;var a=this;this.ajaxRequest.onreadystatechange=function(){a.onXMLLoaded()},this.ajaxRequest.open("GET",this.url,!0),this.ajaxRequest.overrideMimeType&&this.ajaxRequest.overrideMimeType("application/xml"),this.ajaxRequest.send(null)},d.BitmapFontLoader.prototype.onXMLLoaded=function(){if(4===this.ajaxRequest.readyState&&(200===this.ajaxRequest.status||-1===window.location.protocol.indexOf("http"))){var a=this.ajaxRequest.responseXML;if(!a||/MSIE 9/i.test(navigator.userAgent)||navigator.isCocoonJS)if("function"==typeof window.DOMParser){var b=new DOMParser;a=b.parseFromString(this.ajaxRequest.responseText,"text/xml")}else{var c=document.createElement("div");c.innerHTML=this.ajaxRequest.responseText,a=c}var e=this.baseUrl+a.getElementsByTagName("page")[0].getAttribute("file"),f=new d.ImageLoader(e,this.crossorigin);this.texture=f.texture.baseTexture;var g={},h=a.getElementsByTagName("info")[0],i=a.getElementsByTagName("common")[0];g.font=h.getAttribute("face"),g.size=parseInt(h.getAttribute("size"),10),g.lineHeight=parseInt(i.getAttribute("lineHeight"),10),g.chars={};for(var j=a.getElementsByTagName("char"),k=0;k<j.length;k++){var l=parseInt(j[k].getAttribute("id"),10),m=new d.Rectangle(parseInt(j[k].getAttribute("x"),10),parseInt(j[k].getAttribute("y"),10),parseInt(j[k].getAttribute("width"),10),parseInt(j[k].getAttribute("height"),10));g.chars[l]={xOffset:parseInt(j[k].getAttribute("xoffset"),10),yOffset:parseInt(j[k].getAttribute("yoffset"),10),xAdvance:parseInt(j[k].getAttribute("xadvance"),10),kerning:{},texture:d.TextureCache[l]=new d.Texture(this.texture,m)}}var n=a.getElementsByTagName("kerning");for(k=0;k<n.length;k++){var o=parseInt(n[k].getAttribute("first"),10),p=parseInt(n[k].getAttribute("second"),10),q=parseInt(n[k].getAttribute("amount"),10);g.chars[p].kerning[o]=q}d.BitmapText.fonts[g.font]=g;var r=this;f.addEventListener("loaded",function(){r.onLoaded()}),f.load()}},d.BitmapFontLoader.prototype.onLoaded=function(){this.dispatchEvent({type:"loaded",content:this})},d.SpineLoader=function(a,b){d.EventTarget.call(this),this.url=a,this.crossorigin=b,this.loaded=!1},d.SpineLoader.prototype.constructor=d.SpineLoader,d.SpineLoader.prototype.load=function(){var a=this,b=new d.JsonLoader(this.url,this.crossorigin);b.addEventListener("loaded",function(b){a.json=b.content.json,a.onLoaded()}),b.load()},d.SpineLoader.prototype.onLoaded=function(){this.loaded=!0,this.dispatchEvent({type:"loaded",content:this})},d.AbstractFilter=function(a,b){this.passes=[this],this.shaders=[],this.dirty=!0,this.padding=0,this.uniforms=b||{},this.fragmentSrc=a||[]},d.AlphaMaskFilter=function(a){d.AbstractFilter.call(this),this.passes=[this],a.baseTexture._powerOf2=!0,this.uniforms={mask:{type:"sampler2D",value:a},mapDimensions:{type:"2f",value:{x:1,y:5112}},dimensions:{type:"4fv",value:[0,0,0,0]}},a.baseTexture.hasLoaded?(this.uniforms.mask.value.x=a.width,this.uniforms.mask.value.y=a.height):(this.boundLoadedFunction=this.onTextureLoaded.bind(this),a.baseTexture.on("loaded",this.boundLoadedFunction)),this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform sampler2D mask;","uniform sampler2D uSampler;","uniform vec2 offset;","uniform vec4 dimensions;","uniform vec2 mapDimensions;","void main(void) {","   vec2 mapCords = vTextureCoord.xy;","   mapCords += (dimensions.zw + offset)/ dimensions.xy ;","   mapCords.y *= -1.0;","   mapCords.y += 1.0;","   mapCords *= dimensions.xy / mapDimensions;","   vec4 original =  texture2D(uSampler, vTextureCoord);","   float maskAlpha =  texture2D(mask, mapCords).r;","   original *= maskAlpha;","   gl_FragColor =  original;","}"]},d.AlphaMaskFilter.prototype=Object.create(d.AbstractFilter.prototype),d.AlphaMaskFilter.prototype.constructor=d.AlphaMaskFilter,d.AlphaMaskFilter.prototype.onTextureLoaded=function(){this.uniforms.mapDimensions.value.x=this.uniforms.mask.value.width,this.uniforms.mapDimensions.value.y=this.uniforms.mask.value.height,this.uniforms.mask.value.baseTexture.off("loaded",this.boundLoadedFunction)},Object.defineProperty(d.AlphaMaskFilter.prototype,"map",{get:function(){return this.uniforms.mask.value},set:function(a){this.uniforms.mask.value=a}}),d.ColorMatrixFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={matrix:{type:"mat4",value:[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform float invert;","uniform mat4 matrix;","uniform sampler2D uSampler;","void main(void) {","   gl_FragColor = texture2D(uSampler, vTextureCoord) * matrix;","}"]},d.ColorMatrixFilter.prototype=Object.create(d.AbstractFilter.prototype),d.ColorMatrixFilter.prototype.constructor=d.ColorMatrixFilter,Object.defineProperty(d.ColorMatrixFilter.prototype,"matrix",{get:function(){return this.uniforms.matrix.value},set:function(a){this.uniforms.matrix.value=a}}),d.GrayFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={gray:{type:"1f",value:1}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform sampler2D uSampler;","uniform float gray;","void main(void) {","   gl_FragColor = texture2D(uSampler, vTextureCoord);","   gl_FragColor.rgb = mix(gl_FragColor.rgb, vec3(0.2126*gl_FragColor.r + 0.7152*gl_FragColor.g + 0.0722*gl_FragColor.b), gray);","}"]},d.GrayFilter.prototype=Object.create(d.AbstractFilter.prototype),d.GrayFilter.prototype.constructor=d.GrayFilter,Object.defineProperty(d.GrayFilter.prototype,"gray",{get:function(){return this.uniforms.gray.value},set:function(a){this.uniforms.gray.value=a}}),d.DisplacementFilter=function(a){d.AbstractFilter.call(this),this.passes=[this],a.baseTexture._powerOf2=!0,this.uniforms={displacementMap:{type:"sampler2D",value:a},scale:{type:"2f",value:{x:30,y:30}},offset:{type:"2f",value:{x:0,y:0}},mapDimensions:{type:"2f",value:{x:1,y:5112}},dimensions:{type:"4fv",value:[0,0,0,0]}},a.baseTexture.hasLoaded?(this.uniforms.mapDimensions.value.x=a.width,this.uniforms.mapDimensions.value.y=a.height):(this.boundLoadedFunction=this.onTextureLoaded.bind(this),a.baseTexture.on("loaded",this.boundLoadedFunction)),this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform sampler2D displacementMap;","uniform sampler2D uSampler;","uniform vec2 scale;","uniform vec2 offset;","uniform vec4 dimensions;","uniform vec2 mapDimensions;","void main(void) {","   vec2 mapCords = vTextureCoord.xy;","   mapCords += (dimensions.zw + offset)/ dimensions.xy ;","   mapCords.y *= -1.0;","   mapCords.y += 1.0;","   vec2 matSample = texture2D(displacementMap, mapCords).xy;","   matSample -= 0.5;","   matSample *= scale;","   matSample /= mapDimensions;","   gl_FragColor = texture2D(uSampler, vec2(vTextureCoord.x + matSample.x, vTextureCoord.y + matSample.y));","   gl_FragColor.rgb = mix( gl_FragColor.rgb, gl_FragColor.rgb, 1.0);","   vec2 cord = vTextureCoord;","}"]},d.DisplacementFilter.prototype=Object.create(d.AbstractFilter.prototype),d.DisplacementFilter.prototype.constructor=d.DisplacementFilter,d.DisplacementFilter.prototype.onTextureLoaded=function(){this.uniforms.mapDimensions.value.x=this.uniforms.displacementMap.value.width,this.uniforms.mapDimensions.value.y=this.uniforms.displacementMap.value.height,this.uniforms.displacementMap.value.baseTexture.off("loaded",this.boundLoadedFunction)},Object.defineProperty(d.DisplacementFilter.prototype,"map",{get:function(){return this.uniforms.displacementMap.value},set:function(a){this.uniforms.displacementMap.value=a}}),Object.defineProperty(d.DisplacementFilter.prototype,"scale",{get:function(){return this.uniforms.scale.value},set:function(a){this.uniforms.scale.value=a}}),Object.defineProperty(d.DisplacementFilter.prototype,"offset",{get:function(){return this.uniforms.offset.value},set:function(a){this.uniforms.offset.value=a}}),d.PixelateFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={invert:{type:"1f",value:0},dimensions:{type:"4fv",value:new Float32Array([1e4,100,10,10])},pixelSize:{type:"2f",value:{x:10,y:10}}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform vec2 testDim;","uniform vec4 dimensions;","uniform vec2 pixelSize;","uniform sampler2D uSampler;","void main(void) {","   vec2 coord = vTextureCoord;","   vec2 size = dimensions.xy/pixelSize;","   vec2 color = floor( ( vTextureCoord * size ) ) / size + pixelSize/dimensions.xy * 0.5;","   gl_FragColor = texture2D(uSampler, color);","}"]},d.PixelateFilter.prototype=Object.create(d.AbstractFilter.prototype),d.PixelateFilter.prototype.constructor=d.PixelateFilter,Object.defineProperty(d.PixelateFilter.prototype,"size",{get:function(){return this.uniforms.pixelSize.value},set:function(a){this.dirty=!0,this.uniforms.pixelSize.value=a}}),d.BlurXFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={blur:{type:"1f",value:1/512}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform float blur;","uniform sampler2D uSampler;","void main(void) {","   vec4 sum = vec4(0.0);","   sum += texture2D(uSampler, vec2(vTextureCoord.x - 4.0*blur, vTextureCoord.y)) * 0.05;","   sum += texture2D(uSampler, vec2(vTextureCoord.x - 3.0*blur, vTextureCoord.y)) * 0.09;","   sum += texture2D(uSampler, vec2(vTextureCoord.x - 2.0*blur, vTextureCoord.y)) * 0.12;","   sum += texture2D(uSampler, vec2(vTextureCoord.x - blur, vTextureCoord.y)) * 0.15;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y)) * 0.16;","   sum += texture2D(uSampler, vec2(vTextureCoord.x + blur, vTextureCoord.y)) * 0.15;","   sum += texture2D(uSampler, vec2(vTextureCoord.x + 2.0*blur, vTextureCoord.y)) * 0.12;","   sum += texture2D(uSampler, vec2(vTextureCoord.x + 3.0*blur, vTextureCoord.y)) * 0.09;","   sum += texture2D(uSampler, vec2(vTextureCoord.x + 4.0*blur, vTextureCoord.y)) * 0.05;","   gl_FragColor = sum;","}"]},d.BlurXFilter.prototype=Object.create(d.AbstractFilter.prototype),d.BlurXFilter.prototype.constructor=d.BlurXFilter,Object.defineProperty(d.BlurXFilter.prototype,"blur",{get:function(){return this.uniforms.blur.value/(1/7e3)},set:function(a){this.dirty=!0,this.uniforms.blur.value=1/7e3*a}}),d.BlurYFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={blur:{type:"1f",value:1/512}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform float blur;","uniform sampler2D uSampler;","void main(void) {","   vec4 sum = vec4(0.0);","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y - 4.0*blur)) * 0.05;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y - 3.0*blur)) * 0.09;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y - 2.0*blur)) * 0.12;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y - blur)) * 0.15;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y)) * 0.16;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y + blur)) * 0.15;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y + 2.0*blur)) * 0.12;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y + 3.0*blur)) * 0.09;","   sum += texture2D(uSampler, vec2(vTextureCoord.x, vTextureCoord.y + 4.0*blur)) * 0.05;","   gl_FragColor = sum;","}"]},d.BlurYFilter.prototype=Object.create(d.AbstractFilter.prototype),d.BlurYFilter.prototype.constructor=d.BlurYFilter,Object.defineProperty(d.BlurYFilter.prototype,"blur",{get:function(){return this.uniforms.blur.value/(1/7e3)},set:function(a){this.uniforms.blur.value=1/7e3*a}}),d.BlurFilter=function(){this.blurXFilter=new d.BlurXFilter,this.blurYFilter=new d.BlurYFilter,this.passes=[this.blurXFilter,this.blurYFilter]},Object.defineProperty(d.BlurFilter.prototype,"blur",{get:function(){return this.blurXFilter.blur},set:function(a){this.blurXFilter.blur=this.blurYFilter.blur=a}}),Object.defineProperty(d.BlurFilter.prototype,"blurX",{get:function(){return this.blurXFilter.blur},set:function(a){this.blurXFilter.blur=a}}),Object.defineProperty(d.BlurFilter.prototype,"blurY",{get:function(){return this.blurYFilter.blur},set:function(a){this.blurYFilter.blur=a}}),d.InvertFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={invert:{type:"1f",value:1}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform float invert;","uniform sampler2D uSampler;","void main(void) {","   gl_FragColor = texture2D(uSampler, vTextureCoord);","   gl_FragColor.rgb = mix( (vec3(1)-gl_FragColor.rgb) * gl_FragColor.a, gl_FragColor.rgb, 1.0 - invert);","}"]},d.InvertFilter.prototype=Object.create(d.AbstractFilter.prototype),d.InvertFilter.prototype.constructor=d.InvertFilter,Object.defineProperty(d.InvertFilter.prototype,"invert",{get:function(){return this.uniforms.invert.value},set:function(a){this.uniforms.invert.value=a}}),d.SepiaFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={sepia:{type:"1f",value:1}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform float sepia;","uniform sampler2D uSampler;","const mat3 sepiaMatrix = mat3(0.3588, 0.7044, 0.1368, 0.2990, 0.5870, 0.1140, 0.2392, 0.4696, 0.0912);","void main(void) {","   gl_FragColor = texture2D(uSampler, vTextureCoord);","   gl_FragColor.rgb = mix( gl_FragColor.rgb, gl_FragColor.rgb * sepiaMatrix, sepia);","}"]},d.SepiaFilter.prototype=Object.create(d.AbstractFilter.prototype),d.SepiaFilter.prototype.constructor=d.SepiaFilter,Object.defineProperty(d.SepiaFilter.prototype,"sepia",{get:function(){return this.uniforms.sepia.value},set:function(a){this.uniforms.sepia.value=a}}),d.TwistFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={radius:{type:"1f",value:.5},angle:{type:"1f",value:5},offset:{type:"2f",value:{x:.5,y:.5}}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform vec4 dimensions;","uniform sampler2D uSampler;","uniform float radius;","uniform float angle;","uniform vec2 offset;","void main(void) {","   vec2 coord = vTextureCoord - offset;","   float distance = length(coord);","   if (distance < radius) {","       float ratio = (radius - distance) / radius;","       float angleMod = ratio * ratio * angle;","       float s = sin(angleMod);","       float c = cos(angleMod);","       coord = vec2(coord.x * c - coord.y * s, coord.x * s + coord.y * c);","   }","   gl_FragColor = texture2D(uSampler, coord+offset);","}"]},d.TwistFilter.prototype=Object.create(d.AbstractFilter.prototype),d.TwistFilter.prototype.constructor=d.TwistFilter,Object.defineProperty(d.TwistFilter.prototype,"offset",{get:function(){return this.uniforms.offset.value},set:function(a){this.dirty=!0,this.uniforms.offset.value=a}}),Object.defineProperty(d.TwistFilter.prototype,"radius",{get:function(){return this.uniforms.radius.value},set:function(a){this.dirty=!0,this.uniforms.radius.value=a}}),Object.defineProperty(d.TwistFilter.prototype,"angle",{get:function(){return this.uniforms.angle.value},set:function(a){this.dirty=!0,this.uniforms.angle.value=a}}),d.ColorStepFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={step:{type:"1f",value:5}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform sampler2D uSampler;","uniform float step;","void main(void) {","   vec4 color = texture2D(uSampler, vTextureCoord);","   color = floor(color * step) / step;","   gl_FragColor = color;","}"]},d.ColorStepFilter.prototype=Object.create(d.AbstractFilter.prototype),d.ColorStepFilter.prototype.constructor=d.ColorStepFilter,Object.defineProperty(d.ColorStepFilter.prototype,"step",{get:function(){return this.uniforms.step.value},set:function(a){this.uniforms.step.value=a}}),d.DotScreenFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={scale:{type:"1f",value:1},angle:{type:"1f",value:5},dimensions:{type:"4fv",value:[0,0,0,0]}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform vec4 dimensions;","uniform sampler2D uSampler;","uniform float angle;","uniform float scale;","float pattern() {","   float s = sin(angle), c = cos(angle);","   vec2 tex = vTextureCoord * dimensions.xy;","   vec2 point = vec2(","       c * tex.x - s * tex.y,","       s * tex.x + c * tex.y","   ) * scale;","   return (sin(point.x) * sin(point.y)) * 4.0;","}","void main() {","   vec4 color = texture2D(uSampler, vTextureCoord);","   float average = (color.r + color.g + color.b) / 3.0;","   gl_FragColor = vec4(vec3(average * 10.0 - 5.0 + pattern()), color.a);","}"]},d.DotScreenFilter.prototype=Object.create(d.DotScreenFilter.prototype),d.DotScreenFilter.prototype.constructor=d.DotScreenFilter,Object.defineProperty(d.DotScreenFilter.prototype,"scale",{get:function(){return this.uniforms.scale.value},set:function(a){this.dirty=!0,this.uniforms.scale.value=a}}),Object.defineProperty(d.DotScreenFilter.prototype,"angle",{get:function(){return this.uniforms.angle.value},set:function(a){this.dirty=!0,this.uniforms.angle.value=a}}),d.CrossHatchFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={blur:{type:"1f",value:1/512}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform float blur;","uniform sampler2D uSampler;","void main(void) {","    float lum = length(texture2D(uSampler, vTextureCoord.xy).rgb);","    gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);","    if (lum < 1.00) {","        if (mod(gl_FragCoord.x + gl_FragCoord.y, 10.0) == 0.0) {","            gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);","        }","    }","    if (lum < 0.75) {","        if (mod(gl_FragCoord.x - gl_FragCoord.y, 10.0) == 0.0) {","            gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);","        }","    }","    if (lum < 0.50) {","        if (mod(gl_FragCoord.x + gl_FragCoord.y - 5.0, 10.0) == 0.0) {","            gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);","        }","    }","    if (lum < 0.3) {","        if (mod(gl_FragCoord.x - gl_FragCoord.y - 5.0, 10.0) == 0.0) {","            gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);","        }","    }","}"]},d.CrossHatchFilter.prototype=Object.create(d.AbstractFilter.prototype),d.CrossHatchFilter.prototype.constructor=d.BlurYFilter,Object.defineProperty(d.CrossHatchFilter.prototype,"blur",{get:function(){return this.uniforms.blur.value/(1/7e3)},set:function(a){this.uniforms.blur.value=1/7e3*a}}),d.RGBSplitFilter=function(){d.AbstractFilter.call(this),this.passes=[this],this.uniforms={red:{type:"2f",value:{x:20,y:20}},green:{type:"2f",value:{x:-20,y:20}},blue:{type:"2f",value:{x:20,y:-20}},dimensions:{type:"4fv",value:[0,0,0,0]}},this.fragmentSrc=["precision mediump float;","varying vec2 vTextureCoord;","varying vec4 vColor;","uniform vec2 red;","uniform vec2 green;","uniform vec2 blue;","uniform vec4 dimensions;","uniform sampler2D uSampler;","void main(void) {","   gl_FragColor.r = texture2D(uSampler, vTextureCoord + red/dimensions.xy).r;","   gl_FragColor.g = texture2D(uSampler, vTextureCoord + green/dimensions.xy).g;","   gl_FragColor.b = texture2D(uSampler, vTextureCoord + blue/dimensions.xy).b;","   gl_FragColor.a = texture2D(uSampler, vTextureCoord).a;","}"]},d.RGBSplitFilter.prototype=Object.create(d.AbstractFilter.prototype),d.RGBSplitFilter.prototype.constructor=d.RGBSplitFilter,Object.defineProperty(d.RGBSplitFilter.prototype,"angle",{get:function(){return this.uniforms.blur.value/(1/7e3)},set:function(a){this.uniforms.blur.value=1/7e3*a}}),"undefined"!=typeof exports?("undefined"!=typeof module&&module.exports&&(exports=module.exports=d),exports.PIXI=d):"undefined"!=typeof define&&define.amd?define(d):c.PIXI=d}).call(this);
},{}],7:[function(_dereq_,module,exports){
(function (process,Buffer){
/** @license zlib.js 2012 - imaya [ https://github.com/imaya/zlib.js ] The MIT License */(function() {'use strict';function q(b){throw b;}var t=void 0,u=!0;var A="undefined"!==typeof Uint8Array&&"undefined"!==typeof Uint16Array&&"undefined"!==typeof Uint32Array;function E(b,a){this.index="number"===typeof a?a:0;this.m=0;this.buffer=b instanceof(A?Uint8Array:Array)?b:new (A?Uint8Array:Array)(32768);2*this.buffer.length<=this.index&&q(Error("invalid index"));this.buffer.length<=this.index&&this.f()}E.prototype.f=function(){var b=this.buffer,a,c=b.length,d=new (A?Uint8Array:Array)(c<<1);if(A)d.set(b);else for(a=0;a<c;++a)d[a]=b[a];return this.buffer=d};
E.prototype.d=function(b,a,c){var d=this.buffer,f=this.index,e=this.m,g=d[f],k;c&&1<a&&(b=8<a?(G[b&255]<<24|G[b>>>8&255]<<16|G[b>>>16&255]<<8|G[b>>>24&255])>>32-a:G[b]>>8-a);if(8>a+e)g=g<<a|b,e+=a;else for(k=0;k<a;++k)g=g<<1|b>>a-k-1&1,8===++e&&(e=0,d[f++]=G[g],g=0,f===d.length&&(d=this.f()));d[f]=g;this.buffer=d;this.m=e;this.index=f};E.prototype.finish=function(){var b=this.buffer,a=this.index,c;0<this.m&&(b[a]<<=8-this.m,b[a]=G[b[a]],a++);A?c=b.subarray(0,a):(b.length=a,c=b);return c};
var aa=new (A?Uint8Array:Array)(256),J;for(J=0;256>J;++J){for(var N=J,Q=N,ba=7,N=N>>>1;N;N>>>=1)Q<<=1,Q|=N&1,--ba;aa[J]=(Q<<ba&255)>>>0}var G=aa;function R(b,a,c){var d,f="number"===typeof a?a:a=0,e="number"===typeof c?c:b.length;d=-1;for(f=e&7;f--;++a)d=d>>>8^S[(d^b[a])&255];for(f=e>>3;f--;a+=8)d=d>>>8^S[(d^b[a])&255],d=d>>>8^S[(d^b[a+1])&255],d=d>>>8^S[(d^b[a+2])&255],d=d>>>8^S[(d^b[a+3])&255],d=d>>>8^S[(d^b[a+4])&255],d=d>>>8^S[(d^b[a+5])&255],d=d>>>8^S[(d^b[a+6])&255],d=d>>>8^S[(d^b[a+7])&255];return(d^4294967295)>>>0}
var ga=[0,1996959894,3993919788,2567524794,124634137,1886057615,3915621685,2657392035,249268274,2044508324,3772115230,2547177864,162941995,2125561021,3887607047,2428444049,498536548,1789927666,4089016648,2227061214,450548861,1843258603,4107580753,2211677639,325883990,1684777152,4251122042,2321926636,335633487,1661365465,4195302755,2366115317,997073096,1281953886,3579855332,2724688242,1006888145,1258607687,3524101629,2768942443,901097722,1119000684,3686517206,2898065728,853044451,1172266101,3705015759,
2882616665,651767980,1373503546,3369554304,3218104598,565507253,1454621731,3485111705,3099436303,671266974,1594198024,3322730930,2970347812,795835527,1483230225,3244367275,3060149565,1994146192,31158534,2563907772,4023717930,1907459465,112637215,2680153253,3904427059,2013776290,251722036,2517215374,3775830040,2137656763,141376813,2439277719,3865271297,1802195444,476864866,2238001368,4066508878,1812370925,453092731,2181625025,4111451223,1706088902,314042704,2344532202,4240017532,1658658271,366619977,
2362670323,4224994405,1303535960,984961486,2747007092,3569037538,1256170817,1037604311,2765210733,3554079995,1131014506,879679996,2909243462,3663771856,1141124467,855842277,2852801631,3708648649,1342533948,654459306,3188396048,3373015174,1466479909,544179635,3110523913,3462522015,1591671054,702138776,2966460450,3352799412,1504918807,783551873,3082640443,3233442989,3988292384,2596254646,62317068,1957810842,3939845945,2647816111,81470997,1943803523,3814918930,2489596804,225274430,2053790376,3826175755,
2466906013,167816743,2097651377,4027552580,2265490386,503444072,1762050814,4150417245,2154129355,426522225,1852507879,4275313526,2312317920,282753626,1742555852,4189708143,2394877945,397917763,1622183637,3604390888,2714866558,953729732,1340076626,3518719985,2797360999,1068828381,1219638859,3624741850,2936675148,906185462,1090812512,3747672003,2825379669,829329135,1181335161,3412177804,3160834842,628085408,1382605366,3423369109,3138078467,570562233,1426400815,3317316542,2998733608,733239954,1555261956,
3268935591,3050360625,752459403,1541320221,2607071920,3965973030,1969922972,40735498,2617837225,3943577151,1913087877,83908371,2512341634,3803740692,2075208622,213261112,2463272603,3855990285,2094854071,198958881,2262029012,4057260610,1759359992,534414190,2176718541,4139329115,1873836001,414664567,2282248934,4279200368,1711684554,285281116,2405801727,4167216745,1634467795,376229701,2685067896,3608007406,1308918612,956543938,2808555105,3495958263,1231636301,1047427035,2932959818,3654703836,1088359270,
936918E3,2847714899,3736837829,1202900863,817233897,3183342108,3401237130,1404277552,615818150,3134207493,3453421203,1423857449,601450431,3009837614,3294710456,1567103746,711928724,3020668471,3272380065,1510334235,755167117],S=A?new Uint32Array(ga):ga;function ha(){};function ia(b){this.buffer=new (A?Uint16Array:Array)(2*b);this.length=0}ia.prototype.getParent=function(b){return 2*((b-2)/4|0)};ia.prototype.push=function(b,a){var c,d,f=this.buffer,e;c=this.length;f[this.length++]=a;for(f[this.length++]=b;0<c;)if(d=this.getParent(c),f[c]>f[d])e=f[c],f[c]=f[d],f[d]=e,e=f[c+1],f[c+1]=f[d+1],f[d+1]=e,c=d;else break;return this.length};
ia.prototype.pop=function(){var b,a,c=this.buffer,d,f,e;a=c[0];b=c[1];this.length-=2;c[0]=c[this.length];c[1]=c[this.length+1];for(e=0;;){f=2*e+2;if(f>=this.length)break;f+2<this.length&&c[f+2]>c[f]&&(f+=2);if(c[f]>c[e])d=c[e],c[e]=c[f],c[f]=d,d=c[e+1],c[e+1]=c[f+1],c[f+1]=d;else break;e=f}return{index:b,value:a,length:this.length}};function ja(b){var a=b.length,c=0,d=Number.POSITIVE_INFINITY,f,e,g,k,h,l,s,n,m;for(n=0;n<a;++n)b[n]>c&&(c=b[n]),b[n]<d&&(d=b[n]);f=1<<c;e=new (A?Uint32Array:Array)(f);g=1;k=0;for(h=2;g<=c;){for(n=0;n<a;++n)if(b[n]===g){l=0;s=k;for(m=0;m<g;++m)l=l<<1|s&1,s>>=1;for(m=l;m<f;m+=h)e[m]=g<<16|n;++k}++g;k<<=1;h<<=1}return[e,c,d]};function ma(b,a){this.k=na;this.F=0;this.input=A&&b instanceof Array?new Uint8Array(b):b;this.b=0;a&&(a.lazy&&(this.F=a.lazy),"number"===typeof a.compressionType&&(this.k=a.compressionType),a.outputBuffer&&(this.a=A&&a.outputBuffer instanceof Array?new Uint8Array(a.outputBuffer):a.outputBuffer),"number"===typeof a.outputIndex&&(this.b=a.outputIndex));this.a||(this.a=new (A?Uint8Array:Array)(32768))}var na=2,oa={NONE:0,L:1,t:na,X:3},pa=[],T;
for(T=0;288>T;T++)switch(u){case 143>=T:pa.push([T+48,8]);break;case 255>=T:pa.push([T-144+400,9]);break;case 279>=T:pa.push([T-256+0,7]);break;case 287>=T:pa.push([T-280+192,8]);break;default:q("invalid literal: "+T)}
ma.prototype.h=function(){var b,a,c,d,f=this.input;switch(this.k){case 0:c=0;for(d=f.length;c<d;){a=A?f.subarray(c,c+65535):f.slice(c,c+65535);c+=a.length;var e=a,g=c===d,k=t,h=t,l=t,s=t,n=t,m=this.a,p=this.b;if(A){for(m=new Uint8Array(this.a.buffer);m.length<=p+e.length+5;)m=new Uint8Array(m.length<<1);m.set(this.a)}k=g?1:0;m[p++]=k|0;h=e.length;l=~h+65536&65535;m[p++]=h&255;m[p++]=h>>>8&255;m[p++]=l&255;m[p++]=l>>>8&255;if(A)m.set(e,p),p+=e.length,m=m.subarray(0,p);else{s=0;for(n=e.length;s<n;++s)m[p++]=
e[s];m.length=p}this.b=p;this.a=m}break;case 1:var r=new E(A?new Uint8Array(this.a.buffer):this.a,this.b);r.d(1,1,u);r.d(1,2,u);var v=qa(this,f),x,O,y;x=0;for(O=v.length;x<O;x++)if(y=v[x],E.prototype.d.apply(r,pa[y]),256<y)r.d(v[++x],v[++x],u),r.d(v[++x],5),r.d(v[++x],v[++x],u);else if(256===y)break;this.a=r.finish();this.b=this.a.length;break;case na:var D=new E(A?new Uint8Array(this.a.buffer):this.a,this.b),Da,P,U,V,W,qb=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15],ca,Ea,da,Fa,ka,sa=Array(19),
Ga,X,la,B,Ha;Da=na;D.d(1,1,u);D.d(Da,2,u);P=qa(this,f);ca=ra(this.U,15);Ea=ta(ca);da=ra(this.T,7);Fa=ta(da);for(U=286;257<U&&0===ca[U-1];U--);for(V=30;1<V&&0===da[V-1];V--);var Ia=U,Ja=V,I=new (A?Uint32Array:Array)(Ia+Ja),w,K,z,ea,H=new (A?Uint32Array:Array)(316),F,C,L=new (A?Uint8Array:Array)(19);for(w=K=0;w<Ia;w++)I[K++]=ca[w];for(w=0;w<Ja;w++)I[K++]=da[w];if(!A){w=0;for(ea=L.length;w<ea;++w)L[w]=0}w=F=0;for(ea=I.length;w<ea;w+=K){for(K=1;w+K<ea&&I[w+K]===I[w];++K);z=K;if(0===I[w])if(3>z)for(;0<
z--;)H[F++]=0,L[0]++;else for(;0<z;)C=138>z?z:138,C>z-3&&C<z&&(C=z-3),10>=C?(H[F++]=17,H[F++]=C-3,L[17]++):(H[F++]=18,H[F++]=C-11,L[18]++),z-=C;else if(H[F++]=I[w],L[I[w]]++,z--,3>z)for(;0<z--;)H[F++]=I[w],L[I[w]]++;else for(;0<z;)C=6>z?z:6,C>z-3&&C<z&&(C=z-3),H[F++]=16,H[F++]=C-3,L[16]++,z-=C}b=A?H.subarray(0,F):H.slice(0,F);ka=ra(L,7);for(B=0;19>B;B++)sa[B]=ka[qb[B]];for(W=19;4<W&&0===sa[W-1];W--);Ga=ta(ka);D.d(U-257,5,u);D.d(V-1,5,u);D.d(W-4,4,u);for(B=0;B<W;B++)D.d(sa[B],3,u);B=0;for(Ha=b.length;B<
Ha;B++)if(X=b[B],D.d(Ga[X],ka[X],u),16<=X){B++;switch(X){case 16:la=2;break;case 17:la=3;break;case 18:la=7;break;default:q("invalid code: "+X)}D.d(b[B],la,u)}var Ka=[Ea,ca],La=[Fa,da],M,Ma,fa,va,Na,Oa,Pa,Qa;Na=Ka[0];Oa=Ka[1];Pa=La[0];Qa=La[1];M=0;for(Ma=P.length;M<Ma;++M)if(fa=P[M],D.d(Na[fa],Oa[fa],u),256<fa)D.d(P[++M],P[++M],u),va=P[++M],D.d(Pa[va],Qa[va],u),D.d(P[++M],P[++M],u);else if(256===fa)break;this.a=D.finish();this.b=this.a.length;break;default:q("invalid compression type")}return this.a};
function ua(b,a){this.length=b;this.N=a}
var wa=function(){function b(a){switch(u){case 3===a:return[257,a-3,0];case 4===a:return[258,a-4,0];case 5===a:return[259,a-5,0];case 6===a:return[260,a-6,0];case 7===a:return[261,a-7,0];case 8===a:return[262,a-8,0];case 9===a:return[263,a-9,0];case 10===a:return[264,a-10,0];case 12>=a:return[265,a-11,1];case 14>=a:return[266,a-13,1];case 16>=a:return[267,a-15,1];case 18>=a:return[268,a-17,1];case 22>=a:return[269,a-19,2];case 26>=a:return[270,a-23,2];case 30>=a:return[271,a-27,2];case 34>=a:return[272,
a-31,2];case 42>=a:return[273,a-35,3];case 50>=a:return[274,a-43,3];case 58>=a:return[275,a-51,3];case 66>=a:return[276,a-59,3];case 82>=a:return[277,a-67,4];case 98>=a:return[278,a-83,4];case 114>=a:return[279,a-99,4];case 130>=a:return[280,a-115,4];case 162>=a:return[281,a-131,5];case 194>=a:return[282,a-163,5];case 226>=a:return[283,a-195,5];case 257>=a:return[284,a-227,5];case 258===a:return[285,a-258,0];default:q("invalid length: "+a)}}var a=[],c,d;for(c=3;258>=c;c++)d=b(c),a[c]=d[2]<<24|d[1]<<
16|d[0];return a}(),xa=A?new Uint32Array(wa):wa;
function qa(b,a){function c(a,c){var b=a.N,d=[],e=0,f;f=xa[a.length];d[e++]=f&65535;d[e++]=f>>16&255;d[e++]=f>>24;var g;switch(u){case 1===b:g=[0,b-1,0];break;case 2===b:g=[1,b-2,0];break;case 3===b:g=[2,b-3,0];break;case 4===b:g=[3,b-4,0];break;case 6>=b:g=[4,b-5,1];break;case 8>=b:g=[5,b-7,1];break;case 12>=b:g=[6,b-9,2];break;case 16>=b:g=[7,b-13,2];break;case 24>=b:g=[8,b-17,3];break;case 32>=b:g=[9,b-25,3];break;case 48>=b:g=[10,b-33,4];break;case 64>=b:g=[11,b-49,4];break;case 96>=b:g=[12,b-
65,5];break;case 128>=b:g=[13,b-97,5];break;case 192>=b:g=[14,b-129,6];break;case 256>=b:g=[15,b-193,6];break;case 384>=b:g=[16,b-257,7];break;case 512>=b:g=[17,b-385,7];break;case 768>=b:g=[18,b-513,8];break;case 1024>=b:g=[19,b-769,8];break;case 1536>=b:g=[20,b-1025,9];break;case 2048>=b:g=[21,b-1537,9];break;case 3072>=b:g=[22,b-2049,10];break;case 4096>=b:g=[23,b-3073,10];break;case 6144>=b:g=[24,b-4097,11];break;case 8192>=b:g=[25,b-6145,11];break;case 12288>=b:g=[26,b-8193,12];break;case 16384>=
b:g=[27,b-12289,12];break;case 24576>=b:g=[28,b-16385,13];break;case 32768>=b:g=[29,b-24577,13];break;default:q("invalid distance")}f=g;d[e++]=f[0];d[e++]=f[1];d[e++]=f[2];var h,k;h=0;for(k=d.length;h<k;++h)m[p++]=d[h];v[d[0]]++;x[d[3]]++;r=a.length+c-1;n=null}var d,f,e,g,k,h={},l,s,n,m=A?new Uint16Array(2*a.length):[],p=0,r=0,v=new (A?Uint32Array:Array)(286),x=new (A?Uint32Array:Array)(30),O=b.F,y;if(!A){for(e=0;285>=e;)v[e++]=0;for(e=0;29>=e;)x[e++]=0}v[256]=1;d=0;for(f=a.length;d<f;++d){e=k=0;
for(g=3;e<g&&d+e!==f;++e)k=k<<8|a[d+e];h[k]===t&&(h[k]=[]);l=h[k];if(!(0<r--)){for(;0<l.length&&32768<d-l[0];)l.shift();if(d+3>=f){n&&c(n,-1);e=0;for(g=f-d;e<g;++e)y=a[d+e],m[p++]=y,++v[y];break}0<l.length?(s=ya(a,d,l),n?n.length<s.length?(y=a[d-1],m[p++]=y,++v[y],c(s,0)):c(n,-1):s.length<O?n=s:c(s,0)):n?c(n,-1):(y=a[d],m[p++]=y,++v[y])}l.push(d)}m[p++]=256;v[256]++;b.U=v;b.T=x;return A?m.subarray(0,p):m}
function ya(b,a,c){var d,f,e=0,g,k,h,l,s=b.length;k=0;l=c.length;a:for(;k<l;k++){d=c[l-k-1];g=3;if(3<e){for(h=e;3<h;h--)if(b[d+h-1]!==b[a+h-1])continue a;g=e}for(;258>g&&a+g<s&&b[d+g]===b[a+g];)++g;g>e&&(f=d,e=g);if(258===g)break}return new ua(e,a-f)}
function ra(b,a){var c=b.length,d=new ia(572),f=new (A?Uint8Array:Array)(c),e,g,k,h,l;if(!A)for(h=0;h<c;h++)f[h]=0;for(h=0;h<c;++h)0<b[h]&&d.push(h,b[h]);e=Array(d.length/2);g=new (A?Uint32Array:Array)(d.length/2);if(1===e.length)return f[d.pop().index]=1,f;h=0;for(l=d.length/2;h<l;++h)e[h]=d.pop(),g[h]=e[h].value;k=za(g,g.length,a);h=0;for(l=e.length;h<l;++h)f[e[h].index]=k[h];return f}
function za(b,a,c){function d(b){var c=h[b][l[b]];c===a?(d(b+1),d(b+1)):--g[c];++l[b]}var f=new (A?Uint16Array:Array)(c),e=new (A?Uint8Array:Array)(c),g=new (A?Uint8Array:Array)(a),k=Array(c),h=Array(c),l=Array(c),s=(1<<c)-a,n=1<<c-1,m,p,r,v,x;f[c-1]=a;for(p=0;p<c;++p)s<n?e[p]=0:(e[p]=1,s-=n),s<<=1,f[c-2-p]=(f[c-1-p]/2|0)+a;f[0]=e[0];k[0]=Array(f[0]);h[0]=Array(f[0]);for(p=1;p<c;++p)f[p]>2*f[p-1]+e[p]&&(f[p]=2*f[p-1]+e[p]),k[p]=Array(f[p]),h[p]=Array(f[p]);for(m=0;m<a;++m)g[m]=c;for(r=0;r<f[c-1];++r)k[c-
1][r]=b[r],h[c-1][r]=r;for(m=0;m<c;++m)l[m]=0;1===e[c-1]&&(--g[0],++l[c-1]);for(p=c-2;0<=p;--p){v=m=0;x=l[p+1];for(r=0;r<f[p];r++)v=k[p+1][x]+k[p+1][x+1],v>b[m]?(k[p][r]=v,h[p][r]=a,x+=2):(k[p][r]=b[m],h[p][r]=m,++m);l[p]=0;1===e[p]&&d(p)}return g}
function ta(b){var a=new (A?Uint16Array:Array)(b.length),c=[],d=[],f=0,e,g,k,h;e=0;for(g=b.length;e<g;e++)c[b[e]]=(c[b[e]]|0)+1;e=1;for(g=16;e<=g;e++)d[e]=f,f+=c[e]|0,f<<=1;e=0;for(g=b.length;e<g;e++){f=d[b[e]];d[b[e]]+=1;k=a[e]=0;for(h=b[e];k<h;k++)a[e]=a[e]<<1|f&1,f>>>=1}return a};function Aa(b,a){this.input=b;this.b=this.c=0;this.g={};a&&(a.flags&&(this.g=a.flags),"string"===typeof a.filename&&(this.filename=a.filename),"string"===typeof a.comment&&(this.w=a.comment),a.deflateOptions&&(this.l=a.deflateOptions));this.l||(this.l={})}
Aa.prototype.h=function(){var b,a,c,d,f,e,g,k,h=new (A?Uint8Array:Array)(32768),l=0,s=this.input,n=this.c,m=this.filename,p=this.w;h[l++]=31;h[l++]=139;h[l++]=8;b=0;this.g.fname&&(b|=Ba);this.g.fcomment&&(b|=Ca);this.g.fhcrc&&(b|=Ra);h[l++]=b;a=(Date.now?Date.now():+new Date)/1E3|0;h[l++]=a&255;h[l++]=a>>>8&255;h[l++]=a>>>16&255;h[l++]=a>>>24&255;h[l++]=0;h[l++]=Sa;if(this.g.fname!==t){g=0;for(k=m.length;g<k;++g)e=m.charCodeAt(g),255<e&&(h[l++]=e>>>8&255),h[l++]=e&255;h[l++]=0}if(this.g.comment){g=
0;for(k=p.length;g<k;++g)e=p.charCodeAt(g),255<e&&(h[l++]=e>>>8&255),h[l++]=e&255;h[l++]=0}this.g.fhcrc&&(c=R(h,0,l)&65535,h[l++]=c&255,h[l++]=c>>>8&255);this.l.outputBuffer=h;this.l.outputIndex=l;f=new ma(s,this.l);h=f.h();l=f.b;A&&(l+8>h.buffer.byteLength?(this.a=new Uint8Array(l+8),this.a.set(new Uint8Array(h.buffer)),h=this.a):h=new Uint8Array(h.buffer));d=R(s,t,t);h[l++]=d&255;h[l++]=d>>>8&255;h[l++]=d>>>16&255;h[l++]=d>>>24&255;k=s.length;h[l++]=k&255;h[l++]=k>>>8&255;h[l++]=k>>>16&255;h[l++]=
k>>>24&255;this.c=n;A&&l<h.length&&(this.a=h=h.subarray(0,l));return h};var Sa=255,Ra=2,Ba=8,Ca=16;function Y(b,a){this.o=[];this.p=32768;this.e=this.j=this.c=this.s=0;this.input=A?new Uint8Array(b):b;this.u=!1;this.q=Ta;this.K=!1;if(a||!(a={}))a.index&&(this.c=a.index),a.bufferSize&&(this.p=a.bufferSize),a.bufferType&&(this.q=a.bufferType),a.resize&&(this.K=a.resize);switch(this.q){case Ua:this.b=32768;this.a=new (A?Uint8Array:Array)(32768+this.p+258);break;case Ta:this.b=0;this.a=new (A?Uint8Array:Array)(this.p);this.f=this.S;this.z=this.O;this.r=this.Q;break;default:q(Error("invalid inflate mode"))}}
var Ua=0,Ta=1;
Y.prototype.i=function(){for(;!this.u;){var b=Z(this,3);b&1&&(this.u=u);b>>>=1;switch(b){case 0:var a=this.input,c=this.c,d=this.a,f=this.b,e=t,g=t,k=t,h=d.length,l=t;this.e=this.j=0;e=a[c++];e===t&&q(Error("invalid uncompressed block header: LEN (first byte)"));g=e;e=a[c++];e===t&&q(Error("invalid uncompressed block header: LEN (second byte)"));g|=e<<8;e=a[c++];e===t&&q(Error("invalid uncompressed block header: NLEN (first byte)"));k=e;e=a[c++];e===t&&q(Error("invalid uncompressed block header: NLEN (second byte)"));k|=
e<<8;g===~k&&q(Error("invalid uncompressed block header: length verify"));c+g>a.length&&q(Error("input buffer is broken"));switch(this.q){case Ua:for(;f+g>d.length;){l=h-f;g-=l;if(A)d.set(a.subarray(c,c+l),f),f+=l,c+=l;else for(;l--;)d[f++]=a[c++];this.b=f;d=this.f();f=this.b}break;case Ta:for(;f+g>d.length;)d=this.f({B:2});break;default:q(Error("invalid inflate mode"))}if(A)d.set(a.subarray(c,c+g),f),f+=g,c+=g;else for(;g--;)d[f++]=a[c++];this.c=c;this.b=f;this.a=d;break;case 1:this.r(Va,Wa);break;
case 2:Xa(this);break;default:q(Error("unknown BTYPE: "+b))}}return this.z()};
var Ya=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15],Za=A?new Uint16Array(Ya):Ya,$a=[3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258,258,258],ab=A?new Uint16Array($a):$a,bb=[0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0,0,0],cb=A?new Uint8Array(bb):bb,db=[1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577],eb=A?new Uint16Array(db):db,fb=[0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,
10,11,11,12,12,13,13],gb=A?new Uint8Array(fb):fb,hb=new (A?Uint8Array:Array)(288),$,ib;$=0;for(ib=hb.length;$<ib;++$)hb[$]=143>=$?8:255>=$?9:279>=$?7:8;var Va=ja(hb),jb=new (A?Uint8Array:Array)(30),kb,lb;kb=0;for(lb=jb.length;kb<lb;++kb)jb[kb]=5;var Wa=ja(jb);function Z(b,a){for(var c=b.j,d=b.e,f=b.input,e=b.c,g;d<a;)g=f[e++],g===t&&q(Error("input buffer is broken")),c|=g<<d,d+=8;g=c&(1<<a)-1;b.j=c>>>a;b.e=d-a;b.c=e;return g}
function mb(b,a){for(var c=b.j,d=b.e,f=b.input,e=b.c,g=a[0],k=a[1],h,l,s;d<k;){h=f[e++];if(h===t)break;c|=h<<d;d+=8}l=g[c&(1<<k)-1];s=l>>>16;b.j=c>>s;b.e=d-s;b.c=e;return l&65535}
function Xa(b){function a(a,b,c){var d,e,f,g;for(g=0;g<a;)switch(d=mb(this,b),d){case 16:for(f=3+Z(this,2);f--;)c[g++]=e;break;case 17:for(f=3+Z(this,3);f--;)c[g++]=0;e=0;break;case 18:for(f=11+Z(this,7);f--;)c[g++]=0;e=0;break;default:e=c[g++]=d}return c}var c=Z(b,5)+257,d=Z(b,5)+1,f=Z(b,4)+4,e=new (A?Uint8Array:Array)(Za.length),g,k,h,l;for(l=0;l<f;++l)e[Za[l]]=Z(b,3);g=ja(e);k=new (A?Uint8Array:Array)(c);h=new (A?Uint8Array:Array)(d);b.r(ja(a.call(b,c,g,k)),ja(a.call(b,d,g,h)))}
Y.prototype.r=function(b,a){var c=this.a,d=this.b;this.A=b;for(var f=c.length-258,e,g,k,h;256!==(e=mb(this,b));)if(256>e)d>=f&&(this.b=d,c=this.f(),d=this.b),c[d++]=e;else{g=e-257;h=ab[g];0<cb[g]&&(h+=Z(this,cb[g]));e=mb(this,a);k=eb[e];0<gb[e]&&(k+=Z(this,gb[e]));d>=f&&(this.b=d,c=this.f(),d=this.b);for(;h--;)c[d]=c[d++-k]}for(;8<=this.e;)this.e-=8,this.c--;this.b=d};
Y.prototype.Q=function(b,a){var c=this.a,d=this.b;this.A=b;for(var f=c.length,e,g,k,h;256!==(e=mb(this,b));)if(256>e)d>=f&&(c=this.f(),f=c.length),c[d++]=e;else{g=e-257;h=ab[g];0<cb[g]&&(h+=Z(this,cb[g]));e=mb(this,a);k=eb[e];0<gb[e]&&(k+=Z(this,gb[e]));d+h>f&&(c=this.f(),f=c.length);for(;h--;)c[d]=c[d++-k]}for(;8<=this.e;)this.e-=8,this.c--;this.b=d};
Y.prototype.f=function(){var b=new (A?Uint8Array:Array)(this.b-32768),a=this.b-32768,c,d,f=this.a;if(A)b.set(f.subarray(32768,b.length));else{c=0;for(d=b.length;c<d;++c)b[c]=f[c+32768]}this.o.push(b);this.s+=b.length;if(A)f.set(f.subarray(a,a+32768));else for(c=0;32768>c;++c)f[c]=f[a+c];this.b=32768;return f};
Y.prototype.S=function(b){var a,c=this.input.length/this.c+1|0,d,f,e,g=this.input,k=this.a;b&&("number"===typeof b.B&&(c=b.B),"number"===typeof b.M&&(c+=b.M));2>c?(d=(g.length-this.c)/this.A[2],e=258*(d/2)|0,f=e<k.length?k.length+e:k.length<<1):f=k.length*c;A?(a=new Uint8Array(f),a.set(k)):a=k;return this.a=a};
Y.prototype.z=function(){var b=0,a=this.a,c=this.o,d,f=new (A?Uint8Array:Array)(this.s+(this.b-32768)),e,g,k,h;if(0===c.length)return A?this.a.subarray(32768,this.b):this.a.slice(32768,this.b);e=0;for(g=c.length;e<g;++e){d=c[e];k=0;for(h=d.length;k<h;++k)f[b++]=d[k]}e=32768;for(g=this.b;e<g;++e)f[b++]=a[e];this.o=[];return this.buffer=f};
Y.prototype.O=function(){var b,a=this.b;A?this.K?(b=new Uint8Array(a),b.set(this.a.subarray(0,a))):b=this.a.subarray(0,a):(this.a.length>a&&(this.a.length=a),b=this.a);return this.buffer=b};function nb(b){this.input=b;this.c=0;this.G=[];this.R=!1}
nb.prototype.i=function(){for(var b=this.input.length;this.c<b;){var a=new ha,c=t,d=t,f=t,e=t,g=t,k=t,h=t,l=t,s=t,n=this.input,m=this.c;a.C=n[m++];a.D=n[m++];(31!==a.C||139!==a.D)&&q(Error("invalid file signature:"+a.C+","+a.D));a.v=n[m++];switch(a.v){case 8:break;default:q(Error("unknown compression method: "+a.v))}a.n=n[m++];l=n[m++]|n[m++]<<8|n[m++]<<16|n[m++]<<24;a.$=new Date(1E3*l);a.ba=n[m++];a.aa=n[m++];0<(a.n&4)&&(a.W=n[m++]|n[m++]<<8,m+=a.W);if(0<(a.n&Ba)){h=[];for(k=0;0<(g=n[m++]);)h[k++]=
String.fromCharCode(g);a.name=h.join("")}if(0<(a.n&Ca)){h=[];for(k=0;0<(g=n[m++]);)h[k++]=String.fromCharCode(g);a.w=h.join("")}0<(a.n&Ra)&&(a.P=R(n,0,m)&65535,a.P!==(n[m++]|n[m++]<<8)&&q(Error("invalid header crc16")));c=n[n.length-4]|n[n.length-3]<<8|n[n.length-2]<<16|n[n.length-1]<<24;n.length-m-4-4<512*c&&(e=c);d=new Y(n,{index:m,bufferSize:e});a.data=f=d.i();m=d.c;a.Y=s=(n[m++]|n[m++]<<8|n[m++]<<16|n[m++]<<24)>>>0;R(f,t,t)!==s&&q(Error("invalid CRC-32 checksum: 0x"+R(f,t,t).toString(16)+" / 0x"+
s.toString(16)));a.Z=c=(n[m++]|n[m++]<<8|n[m++]<<16|n[m++]<<24)>>>0;(f.length&4294967295)!==c&&q(Error("invalid input size: "+(f.length&4294967295)+" / "+c));this.G.push(a);this.c=m}this.R=u;var p=this.G,r,v,x=0,O=0,y;r=0;for(v=p.length;r<v;++r)O+=p[r].data.length;if(A){y=new Uint8Array(O);for(r=0;r<v;++r)y.set(p[r].data,x),x+=p[r].data.length}else{y=[];for(r=0;r<v;++r)y[r]=p[r].data;y=Array.prototype.concat.apply([],y)}return y};function ob(b){if("string"===typeof b){var a=b.split(""),c,d;c=0;for(d=a.length;c<d;c++)a[c]=(a[c].charCodeAt(0)&255)>>>0;b=a}for(var f=1,e=0,g=b.length,k,h=0;0<g;){k=1024<g?1024:g;g-=k;do f+=b[h++],e+=f;while(--k);f%=65521;e%=65521}return(e<<16|f)>>>0};function pb(b,a){var c,d;this.input=b;this.c=0;if(a||!(a={}))a.index&&(this.c=a.index),a.verify&&(this.V=a.verify);c=b[this.c++];d=b[this.c++];switch(c&15){case rb:this.method=rb;break;default:q(Error("unsupported compression method"))}0!==((c<<8)+d)%31&&q(Error("invalid fcheck flag:"+((c<<8)+d)%31));d&32&&q(Error("fdict flag is not supported"));this.J=new Y(b,{index:this.c,bufferSize:a.bufferSize,bufferType:a.bufferType,resize:a.resize})}
pb.prototype.i=function(){var b=this.input,a,c;a=this.J.i();this.c=this.J.c;this.V&&(c=(b[this.c++]<<24|b[this.c++]<<16|b[this.c++]<<8|b[this.c++])>>>0,c!==ob(a)&&q(Error("invalid adler-32 checksum")));return a};var rb=8;function sb(b,a){this.input=b;this.a=new (A?Uint8Array:Array)(32768);this.k=tb.t;var c={},d;if((a||!(a={}))&&"number"===typeof a.compressionType)this.k=a.compressionType;for(d in a)c[d]=a[d];c.outputBuffer=this.a;this.I=new ma(this.input,c)}var tb=oa;
sb.prototype.h=function(){var b,a,c,d,f,e,g,k=0;g=this.a;b=rb;switch(b){case rb:a=Math.LOG2E*Math.log(32768)-8;break;default:q(Error("invalid compression method"))}c=a<<4|b;g[k++]=c;switch(b){case rb:switch(this.k){case tb.NONE:f=0;break;case tb.L:f=1;break;case tb.t:f=2;break;default:q(Error("unsupported compression type"))}break;default:q(Error("invalid compression method"))}d=f<<6|0;g[k++]=d|31-(256*c+d)%31;e=ob(this.input);this.I.b=k;g=this.I.h();k=g.length;A&&(g=new Uint8Array(g.buffer),g.length<=
k+4&&(this.a=new Uint8Array(g.length+4),this.a.set(g),g=this.a),g=g.subarray(0,k+4));g[k++]=e>>24&255;g[k++]=e>>16&255;g[k++]=e>>8&255;g[k++]=e&255;return g};exports.deflate=ub;exports.deflateSync=vb;exports.inflate=wb;exports.inflateSync=xb;exports.gzip=yb;exports.gzipSync=zb;exports.gunzip=Ab;exports.gunzipSync=Bb;function ub(b,a,c){process.nextTick(function(){var d,f;try{f=vb(b,c)}catch(e){d=e}a(d,f)})}function vb(b,a){var c;c=(new sb(b)).h();a||(a={});return a.H?c:Cb(c)}function wb(b,a,c){process.nextTick(function(){var d,f;try{f=xb(b,c)}catch(e){d=e}a(d,f)})}
function xb(b,a){var c;b.subarray=b.slice;c=(new pb(b)).i();a||(a={});return a.noBuffer?c:Cb(c)}function yb(b,a,c){process.nextTick(function(){var d,f;try{f=zb(b,c)}catch(e){d=e}a(d,f)})}function zb(b,a){var c;b.subarray=b.slice;c=(new Aa(b)).h();a||(a={});return a.H?c:Cb(c)}function Ab(b,a,c){process.nextTick(function(){var d,f;try{f=Bb(b,c)}catch(e){d=e}a(d,f)})}function Bb(b,a){var c;b.subarray=b.slice;c=(new nb(b)).i();a||(a={});return a.H?c:Cb(c)}
function Cb(b){var a=new Buffer(b.length),c,d;c=0;for(d=b.length;c<d;++c)a[c]=b[c];return a};}).call(this); //@ sourceMappingURL=node-zlib.js.map

}).call(this,_dereq_("C:\\Users\\jakec.NEXUSOS\\Documents\\GitHub\\jQuest\\bower_components\\grapefruit\\node_modules\\grunt-browserify\\node_modules\\browserify\\node_modules\\insert-module-globals\\node_modules\\process\\browser.js"),_dereq_("buffer").Buffer)
},{"C:\\Users\\jakec.NEXUSOS\\Documents\\GitHub\\jQuest\\bower_components\\grapefruit\\node_modules\\grunt-browserify\\node_modules\\browserify\\node_modules\\insert-module-globals\\node_modules\\process\\browser.js":5,"buffer":2}],8:[function(_dereq_,module,exports){
var AudioPlayer = _dereq_('./AudioPlayer'),
    inherit = _dereq_('../utils/inherit'),
    support = _dereq_('../utils/support');

//you can only have 1 audio context on a page, so we store one for use in each manager
var __AudioCtx = window.AudioContext || window.webkitAudioContext || window.mozAudioContext,
    __audioctx = support.webAudio ? new __AudioCtx() : null;

/**
 * Grapefruit Audio API, provides an easy interface to use WebAudoiAPI with a fallback to HTML5 Audio
 * The GF Audio API was originally based on [Howler.js](https://github.com/goldfire/howler.js)
 * Generally you will use this via the `game.audio` or `state.audio` properties.
 *
 * @class AudioManager
 * @extends Object
 * @constructor
 * @param game {Game} The game instance this manager belongs to
 * @param parent {AudioManager} The parent audio manager this manager belongs to.
 *      This is used to create a web audio API node heirarchy.
 */
var AudioManager = function(game, parent) {
    /**
     * The game instance this manager belongs to
     *
     * @property game
     * @type Game
     */
    this.game = game;

    /**
     * The parent for this audio manager
     *
     * @property parent
     * @type AudioManager
     */
    this.parent = parent;

    /**
     * Whether the player is muted or not
     *
     * @property muted
     * @type Boolean
     * @default false
     * @private
     */
    this._muted = false;

    /**
     * The master volume of the player
     *
     * @property _volume
     * @type Number
     * @default 1
     * @private
     */
    this._volume = 1;

    /**
     * The Web Audio API context if we are using it
     *
     * @property ctx
     * @type AudioContext
     * @readOnly
     */
    this.ctx = __audioctx;

    /**
     * If we have some way of playing audio
     *
     * @property canPlay
     * @type Boolean
     * @readOnly
     */
    this.canPlay = support.webAudio || support.htmlAudio;

    //if we are using web audio, we need a master gain node
    if(support.webAudio) {
        this.masterGain = this.ctx.createGain ? this.ctx.createGain() : this.ctx.createGainNode();
        this.masterGain.gain.value = 1;
        this.setParent(parent);
    }

    //map of elements to play audio with
    this.sounds = {};
};

inherit(AudioManager, Object, {
    /**
     * Mutes all playing audio
     *
     * @method mute
     * @return {AudioManager} Returns itself.
     * @chainable
     */
    mute: function() {
        return this.setMuted(true);
    },
    /**
     * Unmutes all playing audio
     *
     * @method unmute
     * @return {AudioManager} Returns itself.
     * @chainable
     */
    unmute: function() {
        return this.setMuted(false);
    },
    /**
     * Sets whether or not this manager is muted
     *
     * @method setMuted
     * @return {AudioManager} Returns itself.
     * @chainable
     */
    setMuted: function(m) {
        this._muted = m = !!m;

        if(support.webAudio)
            this.masterGain.gain.value = m ? 0 : this._volume;

        //go through each audio element and mute/unmute them
        for(var key in this.sounds) {
            if(this.sounds.hasOwnProperty(key) && this.sounds[key]._webAudio === false) {
                var player = this.sounds[key];
                //loop through the audio nodes
                for(var i = 0, il = player._nodes.length; i < il; ++i) {
                    player._nodes[i].mute();
                }
            }
        }

        return this;
    },
    /**
     * Sets the parent of this audio manager, if using webAudio this
     * means that we connect to the parent masterGain node and inherit
     * anything that happens to it (such as muting).
     *
     * @method setParent
     * @param parent {AudioManager} The parent to connect to, or `null` to connect to the global context
     * @return {AudioManager} Returns itself.
     * @chainable
     */
    setParent: function(parent) {
        this.parent = parent;

        if(support.webAudio) {
            //attach to parent gain
            if(parent && parent.masterGain) {
                this.masterGain.connect(parent.masterGain);
            }
            //attach to audio context
            else {
                this.masterGain.connect(this.ctx.destination);
            }
        }

        return this;
    },
    /**
     * Attaches an AudioPlayer to this manager, if using webAudio this means
     * that the sound will connect to this masterGain node and inherit anything
     * that happens to it (such as muting).
     *
     * @method attach
     * @param sound {AudioPlayer} The player to attach to this manager
     * @return {AudioPlayer} The newly attached audio player
     */
    attach: function(sound) {
        if(sound._manager !== this) {
            //remove from other manager
            if(sound._manager)
                sound._manager.detach(sound);

            this.sounds[sound.key] = sound;
            sound._manager = this;

            if(support.webAudio) {
                for(var i = 0; i < sound._nodes.length; ++i) {
                    sound._nodes[i].connect(this.masterGain);
                }
            }
        }

        return sound;
    },
    /**
     * Detaches an AudioPlayer from this manager, if using webAudio this means
     * that the sound will disconnect from this masterGain node and stop inheriting
     * anything that happens to it (such as muting).
     *
     * @method detach
     * @param sound {AudioPlayer} The player to detach from this manager
     * @return {AudioPlayer} The detached audio player
     */
    detach: function(sound) {
        if(sound._manager !== this) {
            delete this.sounds[sound.key];
            sound._manager = null;

            if(support.webAudio) {
                for(var i = 0; i < sound._nodes.length; ++i) {
                    sound._nodes[i].disconnect();
                }
            }
        }

        return sound;
    },
    /**
     * Creates a new audio player for a peice of audio
     *
     * @method add
     * @param key {String} The unique cache key for the preloaded audio
     * @param [settings] {Object} All the settings for the audio player
     * @param [settings.volume] {Number} The volume of this audio clip
     * @param [settings.autoplay] {Boolean} Automatically start playing after loaded
     * @param [settings.loop] {Boolean} Replay the audio when it finishes
     * @param [settings.sprite] {Object} A map of string names -> [start, duration] arrays. You can use it to put multiple sounds in one file
     * @param [settings.pos3d] {Array<Number>} 3D coords of where the audio should sound as if it came from (only works with WebAudio)
     * @param [settings.buffer] {Boolean} WebAudio will load the entire file before playing, making this true forces HTML5Audio which will buffer and play
     * @param [settings.format] {String} Force an extension override
     * @return {AudioPlayer} Will return the new audio player, or false if we couldn't determine a compatible url
     */
    add: function(key, settings) {
        //if we can't play audio return false
        if(!this.canPlay) {
            return false;
        }

        var audio = this.game.cache.getAudio(key);

        if(!audio.player)
            audio.player = new AudioPlayer(this, audio, settings);

        return this.sounds[key] = audio.player;
    },
    /**
     * Removes an audio player from the manager
     *
     * @method remove
     * @param key {String} The unique cache key for the preloaded audio
     * @return {AudioPlayer} Will return the audio player removed, or false if none was removed
     */
    remove: function(key) {
        var audio = this.sounds[key];

        if(audio) {
            audio.stop();
        }

        delete this.sounds[key];

        return audio ? audio : false;
    }
});


/**
 * The master volume of all the audio playing
 *
 * @property volume
 * @type Number
 * @default 1
 */
Object.defineProperty(AudioManager.prototype, 'volume', {
    get: function() {
        return this._volume;
    },
    set: function(v) {
        v = parseFloat(v, 10);

        if(!isNaN(v) && v >= 0 && v <= 1) {
            this._volume = v;

            if(support.webAudio)
                this.masterGain.gain.value = v;

            //go through each audio element and change their volume
            for(var key in this.sounds) {
                if(this.sounds.hasOwnProperty(key) && this.sounds[key]._webAudio === false) {
                    var player = this.sounds[key];
                    //loop through the audio nodes
                    for(var i = 0, il = player._nodes.length; i < il; ++i) {
                        player._nodes[i].volume = player._volume * this._volume;
                    }
                }
            }
        }
    }
});

module.exports = AudioManager;

},{"../utils/inherit":69,"../utils/support":70,"./AudioPlayer":9}],9:[function(_dereq_,module,exports){
var AudioPlayer = _dereq_('./AudioPlayer'),
    EventEmitter = _dereq_('../utils/EventEmitter'),
    utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    support = _dereq_('../utils/support');

/**
 * Grapefruit Audio API, provides an easy interface to use HTML5 Audio
 * The GF Audio API was originally based on [Howler.js](https://github.com/goldfire/howler.js)
 *
 * @class AudioPlayer
 * @extends Object
 * @uses EventEmitter
 * @constructor
 * @param manager {AudioManager} AudioManager instance for this audio player
 * @param audio {Object} The preloaded audio file object
 * @param audio.data {ArrayBuffer|Audio} The actual audio data
 * @param audio.webAudio {Boolean} Whether the file is using webAudio or not
 * @param audio.decoded {Boolean} Whether the data has been decoded yet or not
 * @param [settings] {Object} All the settings for this player instance
 * @param [settings.autoplay=false] {Boolean} Whether to automatically start playing the audio file
 * @param [settings.loop=false] {Boolean} Whether the audio should loop or not
 * @param [settings.pos3d] {Array<Number>} The 3d position of the audio to play in the form [x, y, z]
 * @param [settings.sprite] {Object} The audio sprite, if this audio clip has multiple sounds in it.
 *      This object is in the form `{ 'sound': [start, duration] }`, and you can use them with `.play('sound')`.
 */
var AudioPlayer = function(manager, audio, settings) {
    EventEmitter.call(this);

    /**
     * The source of the audio, the actual URL to load up
     *
     * @property src
     * @type String
     */
    this.src = '';

    /**
     * The game instance this player belongs to
     *
     * @property game
     * @type Game
     */
    this.game = manager.game;

    /**
     * The cache key that uniquely identifies this piece of audio
     *
     * @property key
     * @type String
     */
    this.key = audio.key;

    /**
     * Play the audio immediately after loading
     *
     * @property autoplay
     * @type Boolean
     * @default false
     */
    this.autoplay = false;

    /**
     * Override the format determined from the extension with this extension
     *
     * @property format
     * @type String
     * @default null
     */
    this.format = null;

    /**
     * Replay the audio immediately after finishing
     *
     * @property loop
     * @type Boolean
     * @default false
     */
    this.loop = false;

    /**
     * A 3D position where the audio should sound like it is coming from
     *
     * @property pos3d
     * @type Array<Number>
     * @default [0, 0, -0.5]
     */
    this.pos3d = [0, 0, -0.5];

    /**
     * A sound sprite that maps string keys to [start, duration] arrays. These can
     * be used to put multiple sound bits in one file, and play them separately
     *
     * @property sprite
     * @type Object
     * @default {}
     */
    this.sprite = {};

    /**
     * The preloaded audio file object
     *
     * @property _file
     * @type Object
     * @private
     */
    this._file = audio;

    /**
     * The current volume of the player
     *
     * @property _volume
     * @type Number
     * @private
     */
    this._volume = 1;

    /**
     * The full duration of the file to play
     *
     * @property _duration
     * @type Number
     * @private
     */
    this._duration = 0;

    /**
     * Has this player data been loaded?
     *
     * @property _loaded
     * @type Boolean
     * @private
     */
    this._loaded = false;

    /**
     * The manager of the player
     *
     * @property _volum_manager
     * @type AudioManager
     * @private
     */
    this._manager = manager;

    /**
     * Does the browser support WebAudio API
     *
     * @property _webAudio
     * @type Boolean
     * @private
     */
    this._webAudio = support.webAudio;

    /**
     * The actual player nodes, these are either WebAudio Nodes
     * or HTML5 Audio elements.
     *
     * @property _nodes
     * @type Array
     * @private
     */
    this._nodes = [];

    /**
     * Array of timeouts to track end events
     *
     * @property _onendTimer
     * @type Array
     * @private
     */
    this._onendTimer = [];

    //mixin user's settings
    utils.setValues(this, settings);

    if(this._webAudio) {
        this._setupAudioNode();
    }

    this._load();

    /**
     * Fired when the player is ready to play
     *
     * @event ready
     * @param source {String} The source URL that will be used as the audio source
     */

    /**
     * Fired when audio starts playing
     *
     * @event play
     * @param id {String} The id of the node that is used to play the audio
     */

    /**
     * Fired when audio is paused
     *
     * @event paused
     * @param id {String} The id of the node that is paused
     */

    /**
     * Fired when audio finishes playing
     *
     * @event end
     * @param id {String} The id of the node that has finished playing
     */
};

inherit(AudioPlayer, Object, {
    /**
     * Load the audio file for this player, this is called from the ctor
     * there is no reason to call it manually.
     *
     * @method _load
     * @return {AudioPlayer}
     * @private
     */
    _load: function() {
        var self = this,
            audio = this._file;

        //if using web audio, load up the buffer
        if(audio.webAudio) {
            //if not yet decoded, decode before loading the buffer
            if(!audio.decoded) {
                this._manager.ctx.decodeAudioData(audio.data, function(buffer) {
                    if(buffer) {
                        audio.data = buffer;
                        audio.decoded = true;
                        self._loadSoundBuffer(buffer);
                    }
                });
            } else {
                this._loadSoundBuffer(audio.data);
            }
        }
        //otherwise create some Audio nodes
        else {
            //create a new adio node
            var node = audio.data.cloneNode();
            this._nodes.push(node);

            //setup the audio node
            node._pos = 0;
            node.volume = this._manager.muted ? 0 : this._volume * this._manager.volume;

            //setup the event listener to start playing the sound when it has buffered
            this._duration = node.duration;

            //setup a default sprite
            this.sprite._default = [0, node.duration * 1000];

            //check if loaded
            if(!this._loaded) {
                this._loaded = true;
                this.emit('ready', this.src);
            }

            //if autoplay then start it
            if(this.autoplay) {
                this.play();
            }
        }

        return this;
    },
    /**
     * Play a sound from the current time (0 by default).
     *
     * @method play
     * @param [sprite] {String} Plays from the specified position in the sound sprite definition.
     * @param [callback] {Function} Returns the unique playback id for this sound instance.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    play: function(sprite, cb) {
        var self = this;

        if(typeof sprite === 'function') {
            cb = sprite;
            sprite = null;
        }

        //if no sprite specified, use default
        if(!sprite) {
            sprite = '_default';
        }

        //if we haven't loaded yet, wait until we do
        if(!this._loaded) {
            this.on('ready', function() {
                self.play(sprite, cb);
            });

            return this;
        }

        //if the sprite doesn't exist, play nothing
        if(!this.sprite[sprite]) {
            if(typeof cb === 'function') cb();
            return this;
        }

        //get an audio node to use to play
        this._inactiveNode(function(node) {
            var pos = node._pos > 0 ? node._pos : self.sprite[sprite][0] / 1000,
            duration = (self.sprite[sprite][1] / 1000) - node._pos,
            loop = (self.loop || self.sprite[sprite][2]),
            soundId = (typeof cb === 'string') ? cb : (Math.round(Date.now() * Math.random()) + ''),
            timerId;

            node._sprite = sprite;

            //after the audio finishes:
            (function(o) {
                timerId = setTimeout(function() {
                    //if looping restsart it
                    if(!self._webAudio && o.loop) {
                        self.stop(o.id, o.timer).play(o.sprite, o.id);
                    }

                    // set web audio node to paused
                    if(self._webAudio && !o.loop) {
                        self._nodeById(o.id).paused = true;
                    }

                    //end the track if it is HTML audio and a sprite
                    if(!self._webAudio && !o.loop) {
                        self.stop(o.id, o.timer);
                    }

                    //fire off the end event
                    self.emit('end', o.id);
                }, duration * 1000);

                //store the timer
                self._onendTimer.push(timerId);

                //remember which timer to kill
                o.timer = timerId;
            })({
                id: soundId,
                sprite: sprite,
                loop: loop
            });

            //setup webAudio functions
            if(self._webAudio) {
                //set the play id to this node and load into context
                node.id = soundId;
                node.paused = false;
                self._refreshBuffer([loop, pos, duration], soundId);
                self._playStart = self._manager.ctx.currentTime;
                node.gain.value = self._volume;

                if(typeof node.bufferSource.start === 'undefined') {
                    node.bufferSource.noteGrainOn(0, pos, duration);
                } else {
                    node.bufferSource.start(0, pos, duration);
                }
            } else {
                if(node.readyState === 4) {
                    node.id = soundId;
                    node.currentTime = pos;
                    node.muted = self._manager.muted;
                    node.volume = self._volume * self._manager.volume;
                    node.play();
                } else {
                    self._clearEndTimer(timerId);

                    (function() {
                        var sound = self,
                            playSpr = sprite,
                            fn = cb,
                            newNode = node;

                        var evt = function() {
                            sound.play(playSpr, fn);

                            //clear listener
                            newNode.removeEventListener('canplaythrough', evt, false);
                        };
                        newNode.addEventListener('canplaythrough', evt, false);
                    })();

                    return self;
                }
            }

            self.emit('play', soundId);

            if(typeof cb === 'function')
                cb(soundId);
        });

        return this;
    },
    /**
     * Pause playback and save the current position.
     *
     * @method pause
     * @param [id] {String} The play instance ID.
     * @param [timerId] {String} Clear the correct timeout ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    pause: function(id, timerId) {
        var self = this;

        //if we haven't loaded this sound yet, wait until we play it to pause it
        if(!this._loaded) {
            this.on('play', function() {
                self.play(id, timerId);
            });

            return this;
        }

        //clear the onend timer
        this._clearEndTimer(timerId || 0);

        var activeNode = id ? this._nodeById(id) : this._activeNode();
        if(activeNode) {
            if(this._webAudio) {
                //ensure the sound was created
                if(!activeNode.bufferSource)
                    return this;

                activeNode.paused = true;
                activeNode._pos += this._manager.ctx.currentTime - this._playStart;

                if(typeof activeNode.bufferSource.stop === 'undefined') {
                    activeNode.bufferSource.noteOff(0);
                } else {
                    activeNode.bufferSource.stop(0);
                }
            } else {
                activeNode._pos = activeNode.currentTime;
                activeNode.pause();
            }
        }

        this.emit('pause', activeNode ? activeNode.id : id);

        return this;
    },
    /**
     * Stop playback and reset to start.
     *
     * @method stop
     * @param [id] {String} The play instance ID.
     * @param [timerId] {String} Clear the correct timeout ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    stop: function(id, timerId) {
        var self = this;

        //if we haven't loaded this sound yet, wait until we play it to stop it
        if(!this._loaded) {
            this.on('play', function() {
                self.stop(id, timerId);
            });

            return this;
        }

        //clear onend timer
        this._clearEndTimer(timerId || 0);

        var activeNode = id ? this._nodeById(id) : this._activeNode();
        if(activeNode) {
            activeNode._pos = 0;

            if(this._webAudio) {
                if(!activeNode.bufferSource)
                    return this;

                activeNode.paused = true;

                if(typeof activeNode.bufferSource.stop === 'undefined') {
                    activeNode.bufferSource.noteOff(0);
                } else {
                    activeNode.bufferSource.stop(0);
                }
            } else {
                activeNode.pause();
                activeNode.currentTime = 0;
            }
        }

        return this;
    },
    /**
     * Mute this sound.
     *
     * @method mute
     * @param [id] {String} The play instance ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    mute: function(id) {
        return this.setMuted(true, id);
    },
    /**
     * Unmute this sound.
     *
     * @method unmute
     * @param [id] {String} The play instance ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    unmute: function(id) {
        return this.setMuted(false, id);
    },
    /**
     * Set the muted state of this sound.
     *
     * @method setMuted
     * @param muted {Boolean}
     * @param [id] {String} The play instance ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    setMuted: function(muted, id) {
        var self  = this;

        //if we haven't loaded this sound yet, wait until we play it to mute it
        if(!this._loaded) {
            this.on('play', function() {
                self.setMuted(muted, id);
            });

            return this;
        }

        var activeNode = id ? this._nodeById(id) : this._activeNode();
        if(activeNode) {
            if(this._webAudio) {
                activeNode.gain.value = muted ? 0 : this._volume;
            } else {
                activeNode.volume =  muted ? 0 : this._volume;
            }
        }

        return this;
    },
    /**
     * Set the position of playback.
     *
     * @method seek
     * @param pos {Number} The position to move current playback to.
     * @param [id] {String} The play instance ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    seek: function(pos, id) {
        var self = this;

        //if we haven't loaded this sound yet, wait until it is to seek it
        if(!this._loaded) {
            this.on('ready', function() {
                self.seek(pos, id);
            });

            return this;
        }

        //if position is < 0, or invalid, then set to 0
        if(!pos || pos < 0)
            pos = 0;

        var activeNode = id ? this._nodeById(id) : this._activeNode();
        if(activeNode) {
            if(this._webAudio) {
                activeNode._pos = pos;
                this.pause(activeNode.id).play(activeNode._sprite, id);
            } else {
                activeNode.currentTime = pos;
            }
        }

        return this;
    },
    /**
     * Get the position of playback.
     *
     * @method getPosition
     * @param [id] {String} The play instance ID.
     * @return {Number}
     */
    getPosition: function(id) {
        var self = this;

        //if we haven't loaded this sound yet, wait until it is to seek it
        if(!this._loaded) {
            this.on('ready', function() {
                self.getPosition(id);
            });

            return this;
        }

        var activeNode = id ? this._nodeById(id) : this._activeNode();
        if(activeNode) {
            if(this._webAudio) {
                return activeNode._pos + (this._manager.ctx.currentTime - this._playStart);
            } else {
                return activeNode.currentTime;
            }
        }

        return 0;
    },
    /**
     * Fade a currently playing sound between two volumes.
     *
     * @method fade
     * @param from {Number} The volume to fade from (0.0 to 1.0).
     * @param to {Number} The volume to fade to (0.0 to 1.0).
     * @param len {Number} Time in milliseconds to fade.
     * @param [id] {String} The play instance ID.
     * @param [callback] {Function} Fired when the fade is complete.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    fade: function(from, to, len, id, cb) {
        var self = this,
            diff = Math.abs(from - to),
            dir = from > to ? 'dowm' : 'up',
            steps = diff / 0.01,
            stepTime = len / steps;

        if(typeof id === 'function') {
            cb = id;
            id = null;
        }

        //if we haven't loaded this sound yet, wait until it is to seek it
        if(!this._loaded) {
            this.on('ready', function() {
                self.fade(from, to, len, id, cb);
            });

            return this;
        }

        this.setVolume(from, id);

        for(var i = 1; i <= steps; ++i) {
            var change = this._volume + ((dir === 'up' ? 0.01 : -0.01) * i),
                vol = Math.round(1000 * change) / 1000,
                wait = stepTime * i;

            this._doFadeStep(vol, wait, to, id, cb);
        }
    },
    /**
     * Returns the current volume of the player
     *
     * @method getVolume
     * @return {Number} The current volume
     */
    getVolume: function() {
        return this._volume;
    },
    /**
     * Sets the current volume of the player
     *
     * @method setVolume
     * @param vol {Number} The current volume
     * @param [id] {String} The play instance ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    setVolume: function(vol, id) {
        var self = this;

        // make sure volume is a number
        vol = parseFloat(vol);

        //if we haven't loaded this sound yet, wait until we play it to change the volume
        if(!this._loaded) {
            this.on('play', function() {
                self.setVolume(vol, id);
            });

            return this;
        }

        //set the volume
        if(vol >= 0 && vol <= 1) {
            this._volume = vol;

            var activeNode = id ? this._nodeById(id) : this._activeNode();
            if(activeNode) {
                if(this._webAudio) {
                    activeNode.gain.volume = vol;
                } else {
                    activeNode.volume = vol * this._manager.volume;
                }
            }
        }

        return this;
    },
    /**
     * Set the 3D position of the audio source.
     * The most common usage is to set the 'x' position
     * to affect the left/right ear panning. Setting any value higher than
     * 1.0 will begin to decrease the volume of the sound as it moves further away.
     * NOTE: This only works with Web Audio API, HTML5 Audio playback
     * will not be affected.
     *
     * @method setPosition
     * @param x {Number} The x-position of the playback from -1000.0 to 1000.0
     * @param y {Number} The y-position of the playback from -1000.0 to 1000.0
     * @param z {Number} The z-position of the playback from -1000.0 to 1000.0
     * @param [id] {String} The play instance ID.
     * @return {AudioPlayer} Returns itself.
     * @chainable
     */
    setPosition: function(x, y, z, id) {
        var self = this;

        //set a default for the optional 'y' and 'z'
        x = !x ? 0 : x;
        y = !y ? 0 : y;
        z = (!z && z !== 0) ? -0.5 : z;

        //if we haven't loaded this sound yet, wait until we play it to change the position
        if(!this._loaded) {
            this.on('play', function() {
                self.setPosition(x, y, z, id);
            });

            return this;
        }

        if(this._webAudio) {
            var activeNode = id ? this._nodeById(id) : this._activeNode();
            if(activeNode) {
                this.pos3d[0] = x;
                this.pos3d[1] = y;
                this.pos3d[2] = z;
                activeNode.panner.setPosition(x, y, z);
            }
        }

        return this;
    },
    /**
     * Performs a step in the fade transition
     *
     * @method _doFadeStep
     * @private
     */
    _doFadeStep: function(vol, wait, end, id, cb) {
        var self = this;

        setTimeout(function() {
            self.setVolume(vol, id);

            if(vol === end) {
                if(typeof cb === 'function')
                    cb();
            }
        }, wait);
    },
    /**
     * Get an audio node by ID.
     *
     * @method _nodeById
     * @return {AudioPlayer} Audio node.
     * @private
     */
    _nodeById: function(id) {
        var node = this._nodes[0]; //default return value

        //find the node with this ID
        for(var i = 0, il = this._nodes.length; i < il; ++i) {
            if(this._nodes[i].id === id) {
                node = this._nodes[i];
                break;
            }
        }

        return node;
    },
    /**
     * Get the first active audio node.
     *
     * @method _activeNode
     * @return {AudioPlayer} Audio node.
     * @private
     */
    _activeNode: function() {
        var node;

        //find the first playing node
        for(var i = 0, il = this._nodes.length; i < il; ++i) {
            if(!this._nodes[i].paused) {
                node = this._nodes[i];
                break;
            }
        }

        //remove excess inactive nodes
        this._drainPool();

        return node;
    },
    /**
     * Get the first inactive audio node.
     * If there is none, create a new one and add it to the pool.
     *
     * @method _inactiveNode
     * @param cb {Function} callback Function to call when the audio node is ready.
     * @private
     */
    _inactiveNode: function(cb) {
        var node;

        //find first inactive node to recycle
        for(var i = 0, il = this._nodes.length; i < il; ++i) {
            if(this._nodes[i].paused && this._nodes[i].readyState === 4) {
                cb(node = this._nodes[i]);
                break;
            }
        }

        //remove excess inactive nodes
        this._drainPool();

        if(node) return;

        //create new node if there are no inactives
        if(this._webAudio) {
            node = this._setupAudioNode();
            cb(node);
        } else {
            this._load();
            node = this._nodes[this.nodes.length - 1];
            node.addEventListener('loadedmetadata', function() {
                cb(node);
            });
        }
    },
    /**
     * If there are more than 5 inactive audio nodes in the pool, clear out the rest.
     *
     * @method _drainPool
     * @private
     */
    _drainPool: function() {
        var inactive = 0,
            i = 0, il = 0;

        //count inactive nodes
        for(i = 0, il = this._nodes.length; i < il; ++i) {
            if(this._nodes[i].paused) {
                inactive++;
            }
        }

        //remove excess inactive nodes
        for(i = this._nodes.length - 1; i >= 0; --i) {
            if(inactive <= 5)
                break;

            if(this._nodes[i].paused) {
                inactive--;
                this._nodes.splice(i, 1);
            }
        }
    },
    /**
     * Clear 'onend' timeout before it ends.
     *
     * @method _clearEndTimer
     * @param timerId {Number} timerId The ID of the sound to be cancelled.
     * @private
     */
    _clearEndTimer: function(timerId) {
        var timer = this._onendTimer.indexOf(timerId);

        //make sure the timer is cleared
        timer = timer >= 0 ? timer : 0;

        if(this._onendTimer[timer]) {
            clearTimeout(this._onendTimer[timer]);
            this._onendTimer.splice(timer, 1);
        }
    },
    /**
     * Setup the gain node and panner for a Web Audio instance.
     *
     * @method _setupAudioNode
     * @return {Object} The new audio node.
     * @private
     */
    _setupAudioNode: function() {
        var node = this._manager.ctx.createGain ? this._manager.ctx.createGain() : this._manager.ctx.createGainNode();

        this._nodes.push(node);

        //create gain node
        node.gain.value = this._volume;
        node.paused = true;
        node._pos = 0;
        node.readyState = 4;
        node.connect(this._manager.masterGain);

        //create the panner
        node.panner = this._manager.ctx.createPanner();
        node.panner.setPosition(this.pos3d[0], this.pos3d[1], this.pos3d[2]);
        node.panner.connect(node);

        return node;
    },
    /**
     * Finishes loading the Web Audio API sound and fires the loaded event
     *
     * @method loadSound
     * @param buffer {Object} The decoded buffer sound source.
     * @private
     */
    _loadSoundBuffer: function(buffer) {
        this._duration = buffer ? buffer.duration : this._duration;

        //setup a default sprite
        this.sprite._default = [0, this._duration * 1000];

        //fire the load event
        if(!this._loaded) {
            this._loaded = true;
            this.emit('ready', this.src);
        }

        //if autoplay is appropriate
        if(this.autoplay) {
            this.play();
        }
    },
    /**
     * Load the sound back into the buffer source.
     *
     * @method refreshBuffer
     * @param loop {Array} Loop boolean, pos, and duration.
     * @param [id] {String} The play instance ID.
     * @private
     */
    _refreshBuffer: function(loop, id) {
        var node = this._nodeById(id);

        //setup the buffer source for playback
        node.bufferSource = this._manager.ctx.createBufferSource();
        node.bufferSource.buffer = this._file.data;
        node.bufferSource.connect(node.panner);
        node.bufferSource.loop = loop[0];

        if(loop[0]) {
            node.bufferSource.loopStart = loop[1];
            node.bufferSource.loopEnd = loop[1] + loop[2];
        }
    }
});

/**
 * The volume of the audio player
 *
 * @property volume
 * @type Number
 * @default 1
 */
Object.defineProperty(AudioPlayer.prototype, 'volume', {
    get: AudioPlayer.prototype.getVolume,
    set: AudioPlayer.prototype.setVolume
});

module.exports = AudioPlayer;

},{"../utils/EventEmitter":64,"../utils/inherit":69,"../utils/support":70,"../utils/utils":71,"./AudioPlayer":9}],10:[function(_dereq_,module,exports){
var Container = _dereq_('../display/Container'),
    Sprite = _dereq_('../display/Sprite'),
    Rectangle = _dereq_('../geom/Rectangle'),
    Vector = _dereq_('../math/Vector'),
    ObjectPool = _dereq_('../utils/ObjectPool'),
    ObjectFactory = _dereq_('../utils/ObjectFactory'),
    //camera fx
    Close = _dereq_('../fx/camera/Close'),
    Fade = _dereq_('../fx/camera/Fade'),
    Flash = _dereq_('../fx/camera/Flash'),
    Scanlines = _dereq_('../fx/camera/Scanlines'),
    Shake = _dereq_('../fx/camera/Shake'),

    inherit = _dereq_('../utils/inherit'),
    math = _dereq_('../math/math'),
    C = _dereq_('../constants');

/**
 * A basic Camera object that provides some effects. It also will contain the GUI
 * to ensure they are using "screen-coords".
 *
 * @class Camera
 * @extends Container
 * @constructor
 * @param state {State} The game state this camera belongs to
 */
var Camera = function(state) {
    /**
     * The world instance this camera is tied to
     *
     * @property world
     * @type World
     */
    this.world = state.world;

    /**
     * The game instance this camera belongs to
     *
     * @property game
     * @type Game
     */
    this.game = state.game;

    /**
     * The game state this camera belongs to
     *
     * @property state
     * @type State
     */
    this.state = state;

    /**
     * The bounds of that the camera can move to
     *
     * @property bounds
     * @type Rectangle
     * @readOnly
     * @private
     */
    this.bounds = state.world.bounds.clone();

    /**
     * When following a sprite this is the space within the camera that it can move around
     * before the camera moves to track it.
     *
     * @property _deadzone
     * @type Rectangle
     * @readOnly
     * @private
     */
    this._deadzone = null;

    /**
     * The target that the camera will follow
     *
     * @property _target
     * @type Sprite
     * @readOnly
     * @private
     */
    this._target = null;

    /**
     * The target's last position, to cache if we should try to move the camera or not
     *
     * @property _targetPos
     * @type Vector
     * @readOnly
     * @private
     */
    this._targetPos = new Vector();

    /**
     * The size of the camera
     *
     * @property size
     * @type Vector
     * @readOnly
     */
    this.size = new Vector();

    /**
     * Half of the size of the camera
     *
     * @property hSize
     * @type Vector
     * @readOnly
     */
    this.hSize = new Vector();

    /**
     * The container that holds all the GUI items, direct children of Camera are effects
     *
     * @property gui
     * @type Container
     * @readOnly
     */
    this.gui = new Container();

    /**
     * An object factory for creating game objects
     *
     * @property add
     * @type ObjectFactory
     */
    this.add = new ObjectFactory(state, this.gui);

    /**
     * The fxpools for doing camera effects
     *
     * @property fxpools
     * @type Object
     * @private
     * @readOnly
     */
    this.fxpools = {
        flash: new ObjectPool(Flash, this),
        fade: new ObjectPool(Fade, this),
        shake: new ObjectPool(Shake, this),
        scanlines: new ObjectPool(Scanlines, this),
        close: new ObjectPool(Close, this)
    };

    /**
     * Flash the screen with a color. This will cover the screen in a
     * color then fade it out.
     *
     * @method flash
     * @param [color=0xFFFFFF] {Number} The color to flash the screen with
     * @param [duration=1000] {Number} The time it should take (in milliseconds) to fade out
     * @param [alpha=1] {Number} The opacity of the initial flash of color (start opacity)
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Flash} The close effect that was created.
     */

    /**
     * Fade the screen into a color. This will fade into a color that will
     * eventually cover the screen.
     *
     * @method fade
     * @param [color=0xFFFFFF] {Number} The color to fade into
     * @param [duration=1000] {Number} The time it should take (in milliseconds) to fade in
     * @param [alpha=1] {Number} The opacity to fade into (final opacity)
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Fade} The close effect that was created.
     */

    /**
     * Shakes the camera around a bit.
     *
     * @method shake
     * @param [intensity=0.01] {Number} The intensity of the shaking
     * @param [duration=1000] {Number} The amount of time the screen shakes for (in milliseconds)
     * @param [direction=gf.AXIS.BOTH] {gf.AXIS} The axis to shake on
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Shake} The close effect that was created.
     */

    /**
     * Adds arcade-style scanlines to the camera viewport.
     *
     * @method scanlines - color, axis, spacing, thickness, alpha, cb
     * @param [color=0x000000] {Number} The color for the scanlines to be
     * @param [axis=gf.AXIS.HORIZONTAL] {gf.AXIS} The axis to draw the lines on
     * @param [spacing=4] {Number} Number of pixels between each line
     * @param [thickness=1] {Number} Number of pixels thick each line is
     * @param [alpha=0.3] {Number} The opacity of the lines
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Scanlines} The close effect that was created.
     */

    /**
     * Performs a "close" animation that will cover the screen with a color.
     *
     * @method close
     * @param [shape='circle'] {String} The shape to close with, can be either 'ellipse', 'circle', or 'rectangle'
     * @param [duration=1000] {Number} Number of milliseconds for the animation to complete
     * @param [position] {Vector} The position for the animation to close in on, defaults to camera center
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Close} The close effect that was created.
     */

    //Dynamic addition of fx shortcuts
    var self = this;
    Object.keys(this.fxpools).forEach(function(key) {
        self[key] = function() {
            var e = self.fxpools[key].create(),
                args = Array.prototype.slice.call(arguments),
                cb = args.pop();

            if(cb !== undefined && typeof cb !== 'function')
                args.push(cb);

            args.push(this._fxCallback.bind(this, e, key, cb));

            return e.start.apply(e, args);
        };
    });

    Container.call(this);

    //add the gui child
    this.addChild(this.gui);
};

inherit(Camera, Container, {
    /**
     * The base callback for camera FX. This is called at the end of each aniamtion to
     * free the FX class back into the pool.
     *
     * @method _fxCallback
     * @param fx {mixed} The FX instance to free
     * @param type {String} The name of the instance type
     * @param [callback] {Function} The user callback to call.
     * @private
     */
    _fxCallback: function(fx, type, cb) {
        var ret;

        if(typeof cb === 'function')
            ret = cb();

        this.fxpools[type].free(fx);

        return ret;
    },
    /**
     * Follows an sprite with the camera, ensuring they are always center view. You can
     * pass a follow style to change the area an sprite can move around in before we start
     * to move with them.
     *
     * @method follow
     * @param sprite {Sprite} The sprite to follow
     * @param [style=CAMERA_FOLLOW.LOCKON] {CAMERA_FOLLOW} The style of following
     * @return {Camera} Returns itself.
     * @chainable
     */
    follow: function(spr, style) {
        if(!(spr instanceof Sprite))
            return this;

        this._target = spr;
        this._targetPos.set(null, null);

        switch(style) {
            case C.CAMERA_FOLLOW.PLATFORMER:
                var w = this.size.x / 8;
                var h = this.size.y / 3;
                this._deadzone = new Rectangle(
                    (this.size.x - w) / 2,
                    (this.size.y - h) / 2 - (h / 4),
                    w,
                    h
                );
                break;
            case C.CAMERA_FOLLOW.TOPDOWN:
                var sq4 = Math.max(this.size.x, this.size.y) / 4;
                this._deadzone = new Rectangle(
                    (this.size.x - sq4) / 2,
                    (this.size.y - sq4) / 2,
                    sq4,
                    sq4
                );
                break;
            case C.CAMERA_FOLLOW.TOPDOWN_TIGHT:
                var sq8 = Math.max(this.size.x, this.size.y) / 8;
                this._deadzone = new Rectangle(
                    (this.size.x - sq8) / 2,
                    (this.size.y - sq8) / 2,
                    sq8,
                    sq8
                );
                break;
            case C.CAMERA_FOLLOW.LOCKON:
                /* falls through */
            default:
                this._deadzone = null;
                break;
        }

        this.focusSprite(this._target);

        return this;
    },
    /**
     * Stops following any sprites
     *
     * @method unfollow
     * @return {Camera} Returns itself.
     * @chainable
     */
    unfollow: function() {
        this._target = null;
        this._targetPos.set(null, null);
        return this;
    },
    /**
     * Focuses the camera on a sprite.
     *
     * @method focusSprite
     * @param sprite {Sprite} The sprite to focus on
     * @return {Camera} Returns itself.
     * @chainable
     */
    focusSprite: function(spr) {
        var x = spr.position.x,
            y = spr.position.y,
            p = spr.parent;

        //need the transform of the sprite that doesn't take into account
        //the world object. So add up the positions not including the world position.
        while(p && p !== this.world) {
            x += p.position.x;
            y += p.position.y;
            p = p.parent;
        }

        return this.focus(
            //multiple the calculated point by the world scale for this sprite
            x * this.world.scale.x,
            y * this.world.scale.y
        );
    },
    /**
     * Focuses the camera on an x,y position. Ensures that the camera does
     * not go outside the bounds set with setBounds()
     *
     * @method focus
     * @param x {Number|Vector} The x coord to focus on, if a Vector is passed the y param is ignored
     * @param y {Number} The y coord to focus on
     * @return {Camera} Returns itself.
     * @chainable
     */
    focus: function(x, y) {
        y = x.y !== undefined ? x.y : (y || 0);
        x = x.x !== undefined ? x.x : (x || 0);

        //calculate how much we need to pan
        var goToX = x - (this.hSize.x / this.world.worldTransform.a),
            goToY = y - (this.hSize.y / this.world.worldTransform.d),
            dx = goToX + this.world.position.x, //world pos is negative
            dy = goToY + this.world.position.y;

        return this.pan(dx, dy);
    },
    /**
     * Pans the camera around by the x,y amount. Ensures that the camera does
     * not go outside the bounds set with setBounds()
     *
     * @method pan
     * @param x {Number|Vector} The x amount to pan, if a Point is passed the y param is ignored
     * @param y {Number} The y ammount to pan
     * @return {Camera} Returns itself.
     * @chainable
     */
    pan: function(dx, dy) {
        dy = dx.y !== undefined ? dx.y : (dy || 0);
        dx = dx.x !== undefined ? dx.x : (dx || 0);

        if(!dx && !dy) return;

            //world position
        var pos = this.world.position,
            //new world position
            newX = pos.x - dx,
            newY = pos.y - dy,
            b = this.bounds;

        if(b) {
            //check if X movement is illegal
            if(this._outsideBounds(-newX, -pos.y)) {
                dx = (dx < 0 ? b.x : b.right - this.size.x) + pos.x; //how far can we move since dx is too much
            }
            //check if Y movement is illegal
            if(this._outsideBounds(-pos.x, -newY)) {
                dy = (dy < 0 ? b.y : b.bottom - this.size.y) + pos.y;
            }
        }

        if(dx || dy) {
            //prevent NaN
            if(!dx) dx = 0;
            if(!dy) dy = 0;

            this.world.pan(-dx, -dy);
        }

        return this;
    },
    /**
     * Checks if a point is outside the bounds of the camera constraints.
     *
     * @method _outsideBounds
     * @param x {Number} The new X position to test
     * @param y {Number} The new Y position to test
     * @return {Boolean} true if the camera will move outside bounds to go to this point
     * @private
     */
    _outsideBounds: function(x, y) {
        //check if each corner of the camera is within the bounds
        return (
            !this.bounds.contains(x, y) || //top left
            !this.bounds.contains(x, y + this.size.y) || //bottom left
            !this.bounds.contains(x + this.size.x, y) || //top right
            !this.bounds.contains(x + this.size.x, y + this.size.y) //bottom right
        );
    },
    /**
     * Resizes the viewing area, this is called internally by your game instance
     * when you call mygame.resize(). DO NOT CALL THIS DIRECTLY
     *
     * @method resize
     * @private
     * @param w {Number} The new width
     * @param h {Number} The new height
     * @return {Camera} Returns itself.
     * @chainable
     */
    resize: function(w, h) {
        this.size.set(w, h);
        this.hSize.set(
            math.round(this.size.x / 2),
            math.round(this.size.y / 2)
        );

        return this;
    },
    /**
     * Sets the bounds the camera is allowed to go. Usually this is the world's
     * size unless you set it manually.
     *
     * @method constrain
     * @param shape {Rectangle|Polygon|Circle|Ellipse} The shape to constrain the camera into
     * @return {Camera} Returns itself.
     * @chainable
     */
    constrain: function(shape) {
        this.bounds = shape;

        return this;
    },
    /**
     * Removes the constraints of the camera, to allow free movement around the world
     *
     * @method unconstrain
     * @return {Camera} Returns itself.
     * @chainable
     */
    unconstrain: function() {
        this.bounds = null;

        return this;
    },
    /**
     * Called internally every frame. Updates all effects and the follow
     *
     * @method update
     * @param dt {Number} The delta time (in seconds) since the last update
     * @return {Camera} Returns iteself for chainability
     * @private
     */
    update: function(dt) {
        //follow sprite
        if(this._target) {
            var worldTransform = this._target.worldTransform,
                x = worldTransform.tx,
                y = worldTransform.ty;

            if(this._targetPos.x !== x || this._targetPos.y !== y) {
                this._targetPos.set(x, y);

                if(!this._deadzone) {
                    this.focusSprite(this._target);
                } else {
                    var moveX, moveY,
                        dx, dy;

                    moveX = moveY = dx = dy = 0;

                    //check less than
                    dx = x - this._deadzone.x;
                    dy = y - this._deadzone.y;

                    if(dx < 0)
                        moveX = dx;
                    if(dy < 0)
                        moveY = dy;

                    //check greater than
                    dx = x - (this._deadzone.x + this._deadzone.width);
                    dy = y - (this._deadzone.y + this._deadzone.height);

                    if(dx > 0)
                        moveX = dx;
                    if(dy > 0)
                        moveY = dy;

                    this.pan(moveX, moveY);
                }
            }
        }

        //update effects
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var c = this.children[i];
            if(c.update)
                c.update(dt);
        }

        return this;
    }
});

module.exports = Camera;

},{"../constants":11,"../display/Container":16,"../display/Sprite":19,"../fx/camera/Close":23,"../fx/camera/Fade":25,"../fx/camera/Flash":26,"../fx/camera/Scanlines":27,"../fx/camera/Shake":28,"../geom/Rectangle":37,"../math/Vector":47,"../math/math":48,"../utils/ObjectFactory":65,"../utils/ObjectPool":66,"../utils/inherit":69}],11:[function(_dereq_,module,exports){
var constants = {
    /**
     * The types of renderers supported. These are generally passed in to the constructor of
     * a {{#crossLink "Game"}}{{/crossLink}} instance.
     *
     * @class RENDERER
     * @static
     * @final
     */
    RENDERER: {
        /**
         * Represents automatically choosing a renderer
         *
         * @property AUTO
         * @type String
         * @default 'auto'
         * @static
         * @final
         */
        AUTO: 'auto',
        /**
         * Represents the canvas renderer
         *
         * @property CANVAS
         * @type String
         * @default 'canvas'
         * @static
         * @final
         */
        CANVAS: 'canvas',
        /**
         * Represents the webgl renderer
         *
         * @property WEBGL
         * @type String
         * @default 'webgl'
         * @static
         * @final
         */
        WEBGL: 'webgl'
    },

    /**
     * The types of files that the loader supports for types that have multiple formats (like tilemaps).
     *
     * @class FILE_FORMAT
     * @static
     * @final
     */
    FILE_FORMAT: {
        /**
         * Represents the json file type
         *
         * @property JSON
         * @type String
         * @default 'json'
         * @static
         * @final
         */
        JSON: 'json',
        /**
         * Represents the xml file type
         *
         * @property XML
         * @type String
         * @default 'xml'
         * @static
         * @final
         */
        XML: 'xml',
        /**
         * Represents the csv file type
         *
         * @property CSV
         * @type String
         * @default 'csv'
         * @static
         * @final
         */
        CSV: 'csv'
    },

    /**
     * The types of texture atlas file formats that the loader supports.
     *
     * @class ATLAS_FORMAT
     * @static
     * @final
     */
    ATLAS_FORMAT: {
        /**
         * Represents the JSON Array export type of TexturePacker
         *
         * @property JSON_ARRAY
         * @type String
         * @default 'json_array'
         * @static
         * @final
         */
        JSON_ARRAY: 'json_array',
        /**
         * Represents the JSON Hash export type of TexturePacker
         *
         * @property JSON_HASH
         * @type String
         * @default 'json_hash'
         * @static
         * @final
         */
        JSON_HASH: 'json_hash',
        /**
         * Represents the Starling XML format, this export type is supported by TexturePacker
         *
         * @property XML_STARLING
         * @type String
         * @default 'xml_starling'
         * @static
         * @final
         */
        XML_STARLING: 'xml_starling'
    },

    /**
     * The follow types that the camera can execute
     *
     * @class CAMERA_FOLLOW
     * @static
     * @final
     */
    CAMERA_FOLLOW: {
        /**
         * Represents platformer follow style
         *
         * @property PLATFORMER
         * @type Number
         * @default 0
         * @static
         * @final
         */
        PLATFORMER: 0,
        /**
         * Represents topdown follow style
         *
         * @property TOPDOWN
         * @type Number
         * @default 1
         * @static
         * @final
         */
        TOPDOWN: 1,
        /**
         * Represents a tight topdown follow style
         *
         * @property TOPDOWN_TIGHT
         * @type Number
         * @default 2
         * @static
         * @final
         */
        TOPDOWN_TIGHT: 2,
        /**
         * Represents a lockon follow style, this has no deadzone and the camera will
         * follow the target movement exactly.
         *
         * @property LOCKON
         * @type Number
         * @default 3
         * @static
         * @final
         */
        LOCKON: 3
    },

    /**
     * These represent different axis. They are bitwise flags and can be combined together.
     *
     * @class AXIS
     * @static
     * @final
     */
    AXIS: {
        /**
         * Represents no axis, binary value: 0000
         *
         * @property NONE
         * @type Number
         * @default 0
         * @static
         * @final
         */
        NONE: 0,
        /**
         * Represents a horizontal axis, binary value: 0001
         *
         * @property HORIZONTAL
         * @type Number
         * @default 1
         * @static
         * @final
         */
        HORIZONTAL: 1,
        /**
         * Represents a vertical axis, binary value: 0010
         *
         * @property VERTICAL
         * @type Number
         * @default 2
         * @static
         * @final
         */
        VERTICAL: 2,
        /**
         * Represents both axes, binary value: 0011
         *
         * @property VERTICAL
         * @type Number
         * @default 3
         * @static
         * @final
         */
        BOTH: 3
    },

    /**
     * These represent different directions in the world. They are bitwise flags and can be combined together.
     *
     * @class DIRECTION
     * @static
     * @final
     */
    DIRECTION: {
        /**
         * Represents no direction, binary value: 0000
         *
         * @property NONE
         * @type Number
         * @default 0
         * @static
         * @final
         */
        NONE: 0,
        /**
         * Represents left direction, binary value: 0001
         *
         * @property LEFT
         * @type Number
         * @default 1
         * @static
         * @final
         */
        LEFT: 1,
        /**
         * Represents right direction, binary value: 0010
         *
         * @property RIGHT
         * @type Number
         * @default 2
         * @static
         * @final
         */
        RIGHT: 2,
        /**
         * Represents up direction, binary value: 0100
         *
         * @property UP
         * @type Number
         * @default 4
         * @static
         * @final
         */
        UP: 4,
        /**
         * Represents down direction, binary value: 1000
         *
         * @property DOWN
         * @type Number
         * @default 8
         * @static
         * @final
         */
        DOWN: 8,
        /**
         * Represents all directions, binary value: 1111
         *
         * @property ALL
         * @type Number
         * @default 15
         * @static
         * @final
         */
        ALL: 15
    },

    SHAPE: {
        CIRCLE: 1,
        POLYGON: 2,
        RECTANGLE: 3
    }
};

module.exports = constants;

/**
 * The pkg object contains all the grapefruit package information from package.json,
 * inserted at build time.
 *
 * @class pkg
 * @static
 * @final
 */

},{}],12:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit');

/**
 * @class Controls
 * @extends Object
 * @constructor
 * @param game {Game} The game instance this will operate within
 */
var Controls = function(game, actions) {
    this.game = game;
    this.actions = actions || [];
    this.sprites = [];
    this.actionmap = {};

    //create each action object
    for(var i = 0; actions && i < actions.length; ++i) {
        var act = actions[i],
            obj = this.actionmap[act] = {
                callbacks: {}
            };

        //add all the supported bind types (inputs)
        for(var type in Controls.BIND_TYPE) {
            obj[Controls.BIND_TYPE[type]] = []; //binds for an action on an input
            obj.callbacks[type] = null; //the callback for this action on this input
        }
    }
};

inherit(Controls, Object, {
    /**
     * Adds a sprite to this control
     *
     * @method control
     * @param spr {Sprite} The sprite to control
     * @return {Controls} Returns itself
     * @chainable
     */
    control: function(spr) {
        if(this.sprites.indexOf(spr) === -1) {
            this.sprites.push(spr);
        }

        return this;
    },
    /**
     * Removes a sprite from the controls
     *
     * @method release
     * @param spr {Sprite} The sprite to release
     * @return {Controls} Returns itself
     * @chainable
     */
    release: function(spr) {
        var i = this.sprites.indexOf(spr);
        if(i !== -1) {
            this.sprites.splice(i, 1);
        }

        return this;
    },
    setKeys: function(action, keys) {
        this._setBinds(this.game.input.keyboard, Controls.BIND_TYPE.KEYBOARD, action, keys);
    },
    setGpButtons: function(action, buttons) {
        this._setBinds(this.game.input.gamepad.buttons, Controls.BIND_TYPE.GPBUTTON, action, buttons);
    },
    setGpAxis: function(action, axis) {
        this._setBinds(this.game.input.gamepad.sticks, Controls.BIND_TYPE.GPAXIS, action, axis);
    },
    _setBinds: function(input, bindtype, action, keys) {
        var map = this.actionmap[action],
            binds = map[bindtype],
            fn = map.callbacks[bindtype],
            i = 0;

        //remove each bind
        for(i = 0; i < binds.length; ++i) {
            input.off(binds[i], fn);
        }

        //empty binds cache
        binds.length = 0;

        //set new binds
        for(i = 0; i < keys.length; ++i) {
            input.on(keys[i], fn);
            binds.push(keys[i]);
        }
    }
});

Controls.BIND_TYPE = {
    KEYBOARD: 'keyboard',
    GPBUTTON: 'gpbutton',
    GPAXIS: 'gpaxis'
};

module.exports = Controls;

},{"../utils/inherit":69}],13:[function(_dereq_,module,exports){
var Controls = _dereq_('./Controls'),
    Vector = _dereq_('../math/Vector'),
    inherit = _dereq_('../utils/inherit'),
    KEY = _dereq_('../input/Keyboard').KEY,
    BUTTON = _dereq_('../input/gamepad/GamepadButtons').BUTTON,
    AXIS = _dereq_('../input/gamepad/GamepadSticks').AXIS;

/**
 * @class TopDownControls
 * @extends Controls
 * @constructor
 * @param game {Game} The game instance this will operate within
 */
var TopDownControls = function(game, settings) {
    Controls.call(this, game, ['left', 'right', 'up', 'down']);

    //setup callbacks for left
    this.actionmap.left.callbacks[Controls.BIND_TYPE.KEYBOARD] = this.onKey.bind(this, 'left');
    this.actionmap.left.callbacks[Controls.BIND_TYPE.GPBUTTON] = this.onKey.bind(this, 'left');
    this.actionmap.left.callbacks[Controls.BIND_TYPE.GPAXIS] = this.onGpAxis.bind(this, 'left');

    //setup callbacks for right
    this.actionmap.right.callbacks[Controls.BIND_TYPE.KEYBOARD] = this.onKey.bind(this, 'right');
    this.actionmap.right.callbacks[Controls.BIND_TYPE.GPBUTTON] = this.onKey.bind(this, 'right');
    this.actionmap.right.callbacks[Controls.BIND_TYPE.GPAXIS] = this.onGpAxis.bind(this, 'right');

    //setup callbacks for up
    this.actionmap.up.callbacks[Controls.BIND_TYPE.KEYBOARD] = this.onKey.bind(this, 'up');
    this.actionmap.up.callbacks[Controls.BIND_TYPE.GPBUTTON] = this.onKey.bind(this, 'up');
    this.actionmap.up.callbacks[Controls.BIND_TYPE.GPAXIS] = this.onGpAxis.bind(this, 'up');

    //setup callbacks for down
    this.actionmap.down.callbacks[Controls.BIND_TYPE.KEYBOARD] = this.onKey.bind(this, 'down');
    this.actionmap.down.callbacks[Controls.BIND_TYPE.GPBUTTON] = this.onKey.bind(this, 'down');
    this.actionmap.down.callbacks[Controls.BIND_TYPE.GPAXIS] = this.onGpAxis.bind(this, 'down');

    //setup binds
    for(var i = 0; i < this.actions.length; ++i) {
        var act = this.actions[i];

        this.setKeys(act,       settings && settings.keys       ? settings.keys[act]    : TopDownControls.DEFAULT_KEYS[act]);
        this.setGpButtons(act,  settings && settings.buttons    ? settings.buttons[act] : TopDownControls.DEFAULT_BUTTONS[act]);
        this.setGpAxis(act,     settings && settings.axes       ? settings.axes[act]    : TopDownControls.DEFAULT_AXES[act]);
    }

    this.move = {
        left: false,
        right: false,
        up: false,
        down: false,
        lastHorzGpValue: 0,
        lastVertGpValue: 0,
        vec: new Vector(),
        lastVec: new Vector(),
        maxVec: new Vector(1, 1),
        minVec: new Vector(-1, -1),
        dir: {
            left: ['x', -1],
            right: ['x', 1],
            up: ['y', -1],
            down: ['y', 1]
        },
        speed: 100
    };
};

inherit(TopDownControls, Controls, {
    control: function(spr) {
        Controls.prototype.control.call(this, spr);

        spr._phys.system.addControlBody(spr);
        spr._moveVector = new Vector();
    },
    onKey: function(action, evt) {
        if(evt.originalEvent)
            evt.originalEvent.preventDefault();

        // .down is keypressed down
        if(evt.down) {
            if(this.move[action]) return; //skip repeats (holding a key down)

            this.move[action] = true;
            this.move.vec[this.move.dir[action][0]] += this.move.dir[action][1];
        } else {
            this.move[action] = false;
            this.move.vec[this.move.dir[action][0]] -= this.move.dir[action][1];
        }

        this._checkMovement();
    },
    onGpAxis: function(action, evt) {
        if(evt.code === AXIS.LEFT_ANALOGUE_HOR) {
            if(evt.value === 0) {
                if(!this.move.lastHorzGpValue)
                    return;

                this.move.left = false;
                this.move.right = false;
                this.move.vec.x = 0;
            } else if(evt.value > 0) {
                if(this.move.right)
                    return;

                this.move.right = true;
                this.move.vec.x = this.move.dir.right[1];
            } else {
                if(this.move.left)
                    return;

                this.move.left = true;
                this.move.vec.x = this.move.dir.left[1];
            }
            this.move.lastHorzGpValue = evt.value;
        }
        else if(evt.code === AXIS.LEFT_ANALOGUE_VERT) {
            if(evt.value === 0) {
                if(!this.move.lastVertGpValue)
                    return;

                this.move.down = false;
                this.move.up = false;
                this.move.vec.y = 0;
            } else if(evt.value > 0) {
                if(this.move.down)
                    return;

                this.move.down = true;
                this.move.vec.y = this.move.dir.down[1];
            } else {
                if(this.move.up)
                    return;

                this.move.up = true;
                this.move.vec.y = this.move.dir.up[1];
            }
            this.move.lastVertGpValue = evt.value;
        }

        if(!this.move.vec.equals(this.move.lastVec)) {
            this.move.lastVec.copy(this.move.vec);
            this._checkMovement();
        }
    },
    _checkMovement: function() {
        var spr, speed;

        this.move.vec.clamp(this.move.minVec, this.move.maxVec);

        for(var i = 0; i < this.sprites.length; ++i) {
            spr = this.sprites[i];
            speed = spr.moveSpeed || this.move.speed;

            spr._moveVector.x = this.move.vec.x * speed;
            spr._moveVector.y = this.move.vec.y * speed;

            if(spr.locked) return;

            if(spr._setMoveAnimation) spr._setMoveAnimation();

            spr.setVelocity(spr._moveVector);
        }
    }
});

TopDownControls.DEFAULT_KEYS = {
    left:   [KEY.A, KEY.LEFT],
    right:  [KEY.D, KEY.RIGHT],
    up:     [KEY.W, KEY.UP],
    down:   [KEY.S, KEY.DOWN]
};

TopDownControls.DEFAULT_BUTTONS = {
    left:   [BUTTON.PAD_LEFT],
    right:  [BUTTON.PAD_RIGHT],
    up:     [BUTTON.PAD_UP],
    down:   [BUTTON.PAD_DOWN]
};

TopDownControls.DEFAULT_AXES = {
    left:   [AXIS.LEFT_ANALOGUE_HOR, AXIS.RIGHT_ANALOGUE_HOR],
    right:  [AXIS.LEFT_ANALOGUE_HOR, AXIS.RIGHT_ANALOGUE_HOR],
    up:     [AXIS.LEFT_ANALOGUE_VERT, AXIS.RIGHT_ANALOGUE_VERT],
    down:   [AXIS.LEFT_ANALOGUE_VERT, AXIS.RIGHT_ANALOGUE_VERT]
};

module.exports = TopDownControls;

},{"../input/Keyboard":41,"../input/gamepad/GamepadButtons":43,"../input/gamepad/GamepadSticks":44,"../math/Vector":47,"../utils/inherit":69,"./Controls":12}],14:[function(_dereq_,module,exports){
/**
* @license GrapeFruit Game Engine
* Copyright (c) 2012-2014, Chad Engler
*
* GrapeFruit is licensed under the MIT License.
* http://www.opensource.org/licenses/mit-license.php
*
* Known Limiting Features:
*   - Canvas
*       - IE 9+
*       - FF 2+
*       - Chrome 4+
*       - Safari 3.1+
*       - Opera 9+
*
*   - WebGL
*       - IE 11+
*       - FF 4+
*       - Chrome 8+
*       - Safari 6+
*       - Opera 12+
*
*   - Object.create
*       - IE 9+
*       - FF 4+
*       - Chrome 7+
*       - Safari 5+
*       - Opera 12+
*/

var gf = {
    //audio
    AudioManager:       _dereq_('./audio/AudioManager'),
    AudioPlayer:        _dereq_('./audio/AudioPlayer'),

    //camera
    Camera:             _dereq_('./camera/Camera'),

    //controls
    Controls:           _dereq_('./controls/Controls'),
    TopDownControls:    _dereq_('./controls/TopDownControls'),

    //display
    BaseTexture:        _dereq_('./display/BaseTexture'),
    Container:          _dereq_('./display/Container'),
    Graphics:           _dereq_('./display/Graphics'),
    RenderTexture:      _dereq_('./display/RenderTexture'),
    Sprite:             _dereq_('./display/Sprite'),
    SpriteBatch:        _dereq_('./display/SpriteBatch'),
    Texture:            _dereq_('./display/Texture'),
    TilingSprite:       _dereq_('./display/TilingSprite'),

    //fx
    fx: {
        camera: {
            Effect:     _dereq_('./fx/camera/Effect'),
            Close:      _dereq_('./fx/camera/Close'),
            Fade:       _dereq_('./fx/camera/Fade'),
            Flash:      _dereq_('./fx/camera/Flash'),
            Scanlines:  _dereq_('./fx/camera/Scanlines'),
            Shake:      _dereq_('./fx/camera/Shake')
        },
        filters: {
            Filter:     _dereq_('./fx/filters/Filter')
        }
    },

    //game
    Game:               _dereq_('./game/Game'),
    State:              _dereq_('./game/State'),
    StateManager:       _dereq_('./game/StateManager'),
    World:              _dereq_('./game/World'),

    //geometry
    Circle:             _dereq_('./geom/Circle'),
    Ellipse:            _dereq_('./geom/Ellipse'),
    Polygon:            _dereq_('./geom/Polygon'),
    Rectangle:          _dereq_('./geom/Rectangle'),

    //input
    Input:              _dereq_('./input/Input'),
    InputManager:       _dereq_('./input/InputManager'),
    Keyboard:           _dereq_('./input/Keyboard'),
    Gamepad:            _dereq_('./input/Gamepad'),
    GamepadButtons:     _dereq_('./input/gamepad/GamepadButtons'),
    GamepadSticks:      _dereq_('./input/gamepad/GamepadSticks'),
    Pointers:           _dereq_('./input/Pointers'),
    Pointer:            _dereq_('./input/pointer/Pointer'),

    //loader
    Loader:             _dereq_('./loader/Loader'),

    //math
    math:               _dereq_('./math/math'),
    Vector:             _dereq_('./math/Vector'),

    //particles
    ParticleEmitter:    _dereq_('./particles/ParticleEmitter'),
    ParticleSystem:     _dereq_('./particles/ParticleSystem'),

    //physics
    PhysicsSystem:      _dereq_('./physics/PhysicsSystem'),
    PhysicsTarget:      _dereq_('./physics/PhysicsTarget'),

    //text
    BitmapText:         _dereq_('./text/BitmapText'),
    Text:               _dereq_('./text/Text'),

    //tilemap
    Tile:               _dereq_('./tilemap/Tile'),
    Tilelayer:          _dereq_('./tilemap/Tilelayer'),
    Tilemap:            _dereq_('./tilemap/Tilemap'),
    Tileset:            _dereq_('./tilemap/Tileset'),
    ObjectGroup:        _dereq_('./tilemap/ObjectGroup'),

    //utils
    utils:              _dereq_('./utils/utils'),
    support:            _dereq_('./utils/support'),
    inherit:            _dereq_('./utils/inherit'),
    Cache:              _dereq_('./utils/Cache'),
    Clock:              _dereq_('./utils/Clock'),
    Color:              _dereq_('./utils/Color'),
    EventEmitter:       _dereq_('./utils/EventEmitter'),
    ObjectPool:         _dereq_('./utils/ObjectPool'),
    Queue:              _dereq_('./utils/Queue'),
    SpritePool:         _dereq_('./utils/SpritePool'),
    ObjectFactory:      _dereq_('./utils/ObjectFactory'),

    //vendor files
    PIXI:               _dereq_('pixi.js'),
    cp:                 _dereq_('chipmunk')
};

//replace the pixi point with a powerful vector class
gf.PIXI.Point = gf.Vector;

//expose whitelisted pixi filters
var filters = [
    'FilterBlock', 'WebGLFilterManager', 'FilterTexture', 'AlphaMaskFilter',
    'ColorMatrixFilter', 'GrayFilter', 'DisplacementFilter', 'PixelateFilter',
    'BlurXFilter', 'BlurYFilter', 'BlurFilter', 'InvertFilter', 'SepiaFilter',
    'TwistFilter', 'ColorStepFilter', 'DotScreenFilter', 'CrossHatchFilter',
    'RGBSplitFilter'
];

for(var f = 0; f < filters.length; ++f) {
    gf.fx.filters[filters[f]] = gf.PIXI[filters[f]];
}

//copy over constants
var C = _dereq_('./constants');

for(var k in C) {
    gf[k] = C[k];
}

module.exports = gf;

},{"./audio/AudioManager":8,"./audio/AudioPlayer":9,"./camera/Camera":10,"./constants":11,"./controls/Controls":12,"./controls/TopDownControls":13,"./display/BaseTexture":15,"./display/Container":16,"./display/Graphics":17,"./display/RenderTexture":18,"./display/Sprite":19,"./display/SpriteBatch":20,"./display/Texture":21,"./display/TilingSprite":22,"./fx/camera/Close":23,"./fx/camera/Effect":24,"./fx/camera/Fade":25,"./fx/camera/Flash":26,"./fx/camera/Scanlines":27,"./fx/camera/Shake":28,"./fx/filters/Filter":29,"./game/Game":30,"./game/State":31,"./game/StateManager":32,"./game/World":33,"./geom/Circle":34,"./geom/Ellipse":35,"./geom/Polygon":36,"./geom/Rectangle":37,"./input/Gamepad":38,"./input/Input":39,"./input/InputManager":40,"./input/Keyboard":41,"./input/Pointers":42,"./input/gamepad/GamepadButtons":43,"./input/gamepad/GamepadSticks":44,"./input/pointer/Pointer":45,"./loader/Loader":46,"./math/Vector":47,"./math/math":48,"./particles/ParticleEmitter":50,"./particles/ParticleSystem":51,"./physics/PhysicsSystem":52,"./physics/PhysicsTarget":53,"./text/BitmapText":54,"./text/Text":55,"./tilemap/ObjectGroup":56,"./tilemap/Tile":57,"./tilemap/Tilelayer":58,"./tilemap/Tilemap":59,"./tilemap/Tileset":60,"./utils/Cache":61,"./utils/Clock":62,"./utils/Color":63,"./utils/EventEmitter":64,"./utils/ObjectFactory":65,"./utils/ObjectPool":66,"./utils/Queue":67,"./utils/SpritePool":68,"./utils/inherit":69,"./utils/support":70,"./utils/utils":71,"chipmunk":1,"pixi.js":6}],15:[function(_dereq_,module,exports){
/**
 * A texture stores the information that represents an image. All textures have a base texture
 * *This is directly exposing [PIXI.BaseTexture](http://www.goodboydigital.com/pixijs/docs/classes/BaseTexture.html)*
 *
 * @class BaseTexture
 * @uses EventEmitter
 * @constructor
 * @param source {Image|Canvas} the source object (image or canvas)
 */
var BaseTexture = _dereq_('pixi.js').BaseTexture;

module.exports = BaseTexture;

},{"pixi.js":6}],16:[function(_dereq_,module,exports){
var EventEmitter = _dereq_('../utils/EventEmitter'),
    PhysicsTarget = _dereq_('../physics/PhysicsTarget'),
    utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    PIXI = _dereq_('pixi.js');

/**
 * The base display object, that anything being put on the screen inherits from
 * Container or Sprite at some point. This class extends PIXI's DisplayObjectContainer.
 *
 * @class Container
 * @extends [PIXI.DisplayObjectContainer](http://www.goodboydigital.com/pixijs/docs/classes/DisplayObjectContainer.html)
 * @uses EventEmitter
 * @constructor
 */
var Container = function(settings) {
    PIXI.DisplayObjectContainer.call(this);
    EventEmitter.call(this);
    PhysicsTarget.call(this);

    //mixin user's settings
    utils.setValues(this, settings);

    //Add these properties in so that all objects can see them in the docs
    //these properties are inherited from PIXI.DisplayObjectContainer
    //most of these blocks are copied straight from PIXI source

    /**
     * [read-only] The of children of this object.
     * @property children {Array}
     * @readOnly
     */

    /**
     * The coordinate of the object relative to the local coordinates of the parent.
     *
     * @property position
     * @type Point
     */

    /**
     * The scale factor of the object.
     *
     * @property scale
     * @type Point
     */

    /**
     * The rotation of the object in radians.
     *
     * @property rotation
     * @type Number
     */

    /**
     * The opacity of the object.
     *
     * @property alpha
     * @type Number
     */

    /**
     * The visibility of the object.
     *
     * @property visible
     * @type Boolean
     */

    /**
     * [read-only] The display object that contains this display object.
     *
     * @property parent
     * @type DisplayObject
     * @readOnly
     */

    /**
     * [read-only] The stage the display object is connected to, or undefined if it is not connected to the stage.
     *
     * @property stage
     * @type Stage
     * @readOnly
     */

    /**
     * This is the defined area that will pick up mouse / touch events. It is null by default.
     * Setting it is a neat way of optimising the hitTest function that the interactionManager
     * will use (as it will not need to hit test all the children)
     *
     * @property hitArea
     * @type Rectangle|Polygon|Circle|Ellipse
     */

    /**
     * Wether or not the object will handle mouse events
     *
     * @property interactive
     * @type Boolean
     * @default false
     */
};

inherit(Container, PIXI.DisplayObjectContainer, {
    /**
     * Sets the container to visible = true
     *
     * @method show
     * @return {Container} Returns itself.
     * @chainable
     */
    show: function() {
        this.visible = true;
        return this;
    },
    /**
     * Sets the container to visible = false
     *
     * @method hide
     * @return {Container} Returns itself.
     * @chainable
     */
    hide: function() {
        this.visible = false;
        return this;
    },
    /**
     * Adds a child to the container and returns the child
     *
     * @method addChild
     * @param child {Container|Sprite} Any container or sprite
     * @return {Container|Sprite} The child that was added
     */
    addChild: function(child) {
        PIXI.DisplayObjectContainer.prototype.addChild.apply(this, arguments);

        return child;
    },

    /**
     * Adds a child to the object at a specified index. If the index is out of bounds an error will be thrown
     *
     * @method addChildAt
     * @param child {Container|Sprite} Any container or sprite
     * @param index {Number}
     * @return {Container|Sprite} The child that was added
     */
    addChildAt: function(child) {
        PIXI.DisplayObjectContainer.prototype.addChildAt.apply(this, arguments);

        return child;
    },

    /**
     * Removes a child from the object.
     *
     * @method removeChild
     * @param child {Container|Sprite} Any container or sprite
     * @return {Container|Sprite} The child that was added
     */
    removeChild: function(child) {
        PIXI.DisplayObjectContainer.prototype.removeChild.apply(this, arguments);

        return child;
    },

    /**
     * Removes all children from the object.
     *
     * @method removeAllChildren
     * @return {Container} Returns itself.
     * @chainable
     */
    removeAllChildren: function() {
        while(this.children.length) {
            this.removeChild(this.children[0]);
        }

        return this;
    },

    /**
     * Brings a child to the top of the Z pile.
     *
     * @method bringChildToTop
     * @param child {Container|Sprite} Any container or sprite
     * @return {Container|Sprite} The child that was added
     */
    bringChildToTop: function(child) {
        if(child.parent === this) {
            this.addChild(this.removeChild(child));
        }

        return child;
    },

    /**
     * Destroys this object.
     *
     * @method destroy
     */
    destroy: function() {
        this.disablePhysics();
        this.destroyAllChildren();

        if(this.parent)
            this.parent.removeChild(this);
    },

    /**
     * Destroys all the children of the object.
     *
     * @method destroyAllChildren
     * @return {Container} Returns itself.
     * @chainable
     */
    destroyAllChildren: function() {
        while(this.children.length) {
            if(this.children[0].destroy) {
                this.children[0].destroy();
            } else {
                this.removeChild(this.children[0]);
            }
        }

        return this;
    }
});

module.exports = Container;

},{"../physics/PhysicsTarget":53,"../utils/EventEmitter":64,"../utils/inherit":69,"../utils/utils":71,"pixi.js":6}],17:[function(_dereq_,module,exports){
/**
 * The Graphics class contains a set of methods that you can use to create primitive shapes and lines.
 * It is important to know that with the webGL renderer only simple polys can be filled at this stage
 * Complex polys will not be filled. Heres an example of a
 * [complex polygon](http://www.goodboydigital.com/wp-content/uploads/2013/06/complexPolygon.png).
 * *This is directly exposing [PIXI.Graphics](http://www.goodboydigital.com/pixijs/docs/classes/Graphics.html)*
 *
 * @class Graphics
 * @extends [PIXI.DisplayObjectContainer](http://www.goodboydigital.com/pixijs/docs/classes/DisplayObjectContainer.html)
 * @constructor
 */
var Graphics = _dereq_('pixi.js').Graphics;

module.exports = Graphics;

},{"pixi.js":6}],18:[function(_dereq_,module,exports){
/**
 * A RenderTexture is a special texture that allows any pixi displayObject to be rendered to it.
 *
 * __Hint__: All DisplayObjects (exmpl. Sprites) that renders on RenderTexture should be preloaded.
 * Otherwise black rectangles will be drawn instead.
 *
 * RenderTexture takes snapshot of DisplayObject passed to render method. If DisplayObject is passed to render method, position and rotation of it will be ignored. For example:
 *
 * ```
 * var renderTexture = new gf.RenderTexture(800, 600);
 * var sprite = gf.Sprite(texture);
 * sprite.position.x = 800/2;
 * sprite.position.y = 600/2;
 * sprite.anchor.x = 0.5;
 * sprite.anchor.y = 0.5;
 * renderTexture.render(sprite);
 * ```
 *
 * Sprite in this case will be rendered to 0,0 position. To render this sprite at center Container should be used:
 *
 * ```
 * var doc = new gf.Container();
 * doc.addChild(sprite);
 * renderTexture.render(doc);  // Renders to center of renderTexture
 * ```
 *
 * *This is directly exposing [PIXI.RenderTexture](http://www.goodboydigital.com/pixijs/docs/classes/RenderTexture.html)*
 *
 * @class RenderTexture
 * @extends Texture
 * @constructor
 * @param width {Number} The width of the render texture
 * @param height {Number} The height of the render texture
 */
var RenderTexture = _dereq_('pixi.js').RenderTexture;

module.exports = RenderTexture;

},{"pixi.js":6}],19:[function(_dereq_,module,exports){
var EventEmitter = _dereq_('../utils/EventEmitter'),
    Rectangle = _dereq_('../geom/Rectangle'),
    PhysicsTarget = _dereq_('../physics/PhysicsTarget'),
    Clock = _dereq_('../utils/Clock'),
    inherit = _dereq_('../utils/inherit'),
    Texture = _dereq_('./Texture'),
    utils = _dereq_('../utils/utils'),
    PIXI = _dereq_('pixi.js');

/**
 * The base Sprite class. This class is the base for all images on the screen. This class extends PIXI's Sprite.
 *
 * @class Sprite
 * @extends [PIXI.Sprite](http://www.goodboydigital.com/pixijs/docs/classes/Sprite.html)
 * @uses EventEmitter
 * @uses PhysicsTarget
 * @constructor
 * @param textures {Texture|Array<Texture>|Object} The texture for the sprite to display, an array of texture to animation through, or an animation object.
 *      The later looks like: `{ animationName: { frames: [frame1, frame2], frameTime: 0.5, loop: false } }` where each frame is a Texture object
 * @param [frameTime] {Number} The frame time of the animations (can be overriden on a specific animations)
 * @param [start] {String} The animation to start with, defaults to the first found key otherwise
 * @example
 *      var spr = new gf.Sprite(texture);
 */
var Sprite = function(anims, frameTime, start) {
    EventEmitter.call(this);
    PhysicsTarget.call(this);

    if(!anims) {
        anims = Texture.__default;
    }

    //parse tx into correct format
    if(anims instanceof Texture) {
        anims = { _default: { frames: [anims] } };
        start = '_default';
    }
    else if(anims instanceof Array) {
        anims = { _default: { frames: anims } };
        start = '_default';
    } else {
        //massage animations into full format
        for(var a in anims) {
            if(start === undefined)
                start = a;

            var anim = anims[a];

            if(anim instanceof Array)
                anims[a] = { name: a, frames: anim };
            else if(anim instanceof Texture)
                anims[a] = { name: a, frames: [anim] };
            else
                anims[a].name = a;
        }
    }

    PIXI.Sprite.call(this, anims[start].frames[0]);

    /**
     * The name of the sprite
     *
     * @property name
     * @type String
     * @default ''
     */
    this.name = '';

    /**
     * The lifetime of the sprite. Once it reaches 0 (after being set)
     * the sprite's visible property is set to false, so that it will
     * no longer be rendered. NOT YET IMPLEMENTED
     *
     * @property lifetime
     * @type Number
     * @default Infinity
     * @private
     */
    this.lifespan = Infinity;

    /**
     * The animation frameTime for this sprite
     *
     * @property frameTime
     * @type Number
     * @default 1
     */
    this.frameTime = frameTime || 250;

    /**
     * Whether or not to loop the animations. This can be overriden
     * on a per-animation level
     *
     * @property loop
     * @type Boolean
     * @default false
     */
    this.loop = false;

    /**
     * The registerd animations for this AnimatedSprite
     *
     * @property animations
     * @type Object
     * @readOnly
     */
    this.animations = anims;

    /**
     * The currently playing animation
     *
     * @property currentAnimation
     * @type String
     * @readOnly
     */
    this.currentAnimation = start;

    /**
     * The current frame being shown
     *
     * @property currentFrame
     * @type Number
     * @readOnly
     */
    this.currentFrame = 0;

    /**
     * Whether or not the animation is currently playing
     *
     * @property playing
     * @type Boolean
     * @readOnly
     */
    this.playing = false;

    this.hitArea = this.hitArea || new Rectangle(0, 0, this.width, this.height);

    //show first frame
    this.goto(0, this.currentAnimation);

    /**
     * Fired when a new frame of the running animation is shown
     *
     * @event frame
     * @param animation {String} The animation name that is playing
     * @param frameId {Number} The frame that is being shown
     */

    /**
     * Fired when the running animation completes
     *
     * @event complete
     * @param animation {String} The animation that has completed
     */

    this._clock = new Clock();
    this._clock.start();
    this._frameTime = 0;
};

inherit(Sprite, PIXI.Sprite, {
    enableInteractivity: function(game) {
        game.input.pointers.watchSprite(this);
        return this;
    },
    /**
     * Sets the sprite to visible = true
     *
     * @method show
     * @return {Sprite} Returns itself.
     * @chainable
     */
    show: function() {
        this.visible = true;
        return this;
    },
    /**
     * Sets the sprite to visible = false
     *
     * @method hide
     * @return {Sprite} Returns itself.
     * @chainable
     */
    hide: function() {
        this.visible = false;
        return this;
    },
    /**
     * Creates a new Sprite instance with the same values as this one
     *
     * @method clone
     * @return {Sprite} Returns the new sprite
     */
    clone: function(spr) {
        //make a copy of our animations object
        var anims = utils.extend(true, {}, this.animations);

        spr = spr || new Sprite(anims, this.frameTime, this.currentAnimation);

        spr.name = this.name;
        spr.loop = this.loop;
        spr.currentFrame = this.currentFrame;
        spr.playing = this.playing;
        spr.hitArea = this.hitArea.clone();

        spr.blendMode = this.blendMode;

        spr.anchor.x = this.anchor.x;
        spr.anchor.y = this.anchor.y;

        spr.position.x = this.position.x;
        spr.position.y = this.position.y;

        spr.scale.x = this.scale.x;
        spr.scale.y = this.scale.y;

        spr.pivot.x = this.pivot.x;
        spr.pivot.y = this.pivot.y;

        spr.rotation = this.rotation;
        spr.alpha = this.alpha;
        spr.visible = this.visible;
        spr.buttonMode = this.buttonMode;
        spr.renderable = this.renderable;

        //Don't copy some pixi stuff

        //spr.children = this.children;
        //spr.parent = this.parent;
        //spr.stage = this.stage;
        //spr.worldAlpha = this.worldAlpha;
        //spr._interactive = this._interactive;
        //spr.worldTransform = this.worldTransform;
        //spr.localTransform = this.localTransform;
        //spr.color = this.color;
        //spr.dynamic = this.dynamic;

        return spr;
    },
    /**
     * Adds a new animation to this animated sprite
     *
     * @method addAnimation
     * @param name {String} The string name of the animation
     * @param frames {Array<Texture>} The array of texture frames
     * @param [frameTime] {Number} The animation frame time
     * @param [loop] {Boolean} Loop the animation or not
     * @return {Sprite} Returns itself.
     * @chainable
     */
    addAnimation: function(name, frames, frameTime, loop) {
        if(typeof name === 'object') {
            this.animations[name.name] = name;
        } else {
            this.animations[name] = {
                name: name,
                frames: frames,
                frameTime: frameTime,
                loop: loop
            };
        }

        return this;
    },
    /**
     * Goes to a frame and starts playing the animation from there. You can optionally
     * pass the name of a new aniamtion to start playing.
     *
     * @method goto
     * @param frame {Number} The index of the frame to start on
     * @param [name] {String} The string name of the animation to go to
     * @return {Sprite} Returns itself.
     * @chainable
     */
    goto: function(frame, anim) {
        if(typeof frame === 'string') {
            anim = frame;
            frame = 0;
        }

        this.currentFrame = frame || 0;

        if(anim) {
            this.currentAnimation = anim;
        }

        this.setTexture(this.animations[this.currentAnimation].frames[this.currentFrame]);
        this.emit('frame', this.currentAnimation, this.currentFrame);

        return this;
    },
    /**
     * Starts playing the currently active animation
     *
     * @method play
     * @return {Sprite} Returns itself.
     * @chainable
     */
    play: function() {
        this.playing = true;
        return this;
    },
    /**
     * Stops playing the currently active animation
     *
     * @method stop
     * @return {Sprite} Returns itself.
     * @chainable
     */
    stop: function() {
        this.playing = false;
        return this;
    },
    /**
     * Removes this sprite from the stage and the physics system
     *
     * @method destroy
     */
    destroy: function() {
        this.stop();
        this.disablePhysics();

        if(this.parent)
            this.parent.removeChild(this);

        this.name = null;
        this.lifetime = null;
        this.frameTime = null;
        this.loop = null;
        this.animations = null;
        this.currentAnimation = null;
        this.currentFrame = null;
        this.playing = null;
        this.hitArea = null;
    },
    /**
     * Called by PIXI to update our textures and do the actual animation
     *
     * @method updateTransform
     * @private
     */
    updateTransform: function() {
        PIXI.Sprite.prototype.updateTransform.call(this);

        //get delta
        var dt = this._clock.getDelta();

        //if not playing, quit
        if(!this.playing) return;

        //get anim and test time
        var anim = this.animations[this.currentAnimation];

        this._frameTime += dt * 1000;

        if(this._frameTime < (anim.frameTime || 500)) return;

        //update animation to next frame
        var loop = anim.loop !== undefined ? anim.loop : this.loop,
            frame = ++this.currentFrame;

        this._frameTime = 0;

        if(frame < anim.frames.length) {
            this.setTexture(anim.frames[frame]);
            this.emit('frame', this.currentAnimation, frame);
        }
        else {
            if(loop) {
                this.goto(0);
            } else {
                this.stop();
                this.emit('complete', this.currentAnimation);
            }
        }
    }
});

module.exports = Sprite;

//Add event echos
/*
['click', 'mousedown', 'mouseup', 'mouseupoutside', 'mouseover', 'mouseout', 'mousemove', 'tap', 'touchstart', 'touchend', 'touchendoutside'].forEach(function(evtname) {
    Sprite.prototype[evtname] = module.exports = function(e) {
        this.emit(evtname, e);
    };
});
*/

/*
 * MOUSE Callbacks
 */

/**
 * A callback that is used when the users clicks on the sprite with their mouse
 *
 * @event click
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the user clicks the mouse down over the sprite
 *
 * @event mousedown
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the user releases the mouse that was over the sprite
 * for this callback to be fired the mouse must have been pressed down over the sprite
 *
 * @event mouseup
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the user releases the mouse that was over the sprite but is no longer over the sprite
 * for this callback to be fired, The touch must have started over the sprite
 *
 * @event mouseupoutside
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the users mouse rolls over the sprite
 *
 * @event mouseover
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the users mouse leaves the sprite
 *
 * @event mouseout
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the user moves the mouse while over the sprite
 *
 * @event mousemove
 * @param interactionData {InteractionData}
 */

/*
 * TOUCH Callbacks
 */

/**
 * A callback that is used when the users taps on the sprite with their finger
 * basically a touch version of click
 *
 * @event tap
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the user touch's over the sprite
 *
 * @event touchstart
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the user releases a touch over the sprite
 *
 * @event touchend
 * @param interactionData {InteractionData}
 */

/**
 * A callback that is used when the user releases the touch that was over the sprite
 * for this callback to be fired, The touch must have started over the sprite
 *
 * @event touchendoutside
 * @param interactionData {InteractionData}
 */

},{"../geom/Rectangle":37,"../physics/PhysicsTarget":53,"../utils/Clock":62,"../utils/EventEmitter":64,"../utils/inherit":69,"../utils/utils":71,"./Texture":21,"pixi.js":6}],20:[function(_dereq_,module,exports){
var EventEmitter = _dereq_('../utils/EventEmitter'),
    utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    PIXI = _dereq_('pixi.js');

/**
 * The same as a container, but offers improvements in render speed for single-level
 * deep sprite children. Useful for containers like tilemap layers, bitmap text container,
 * particle containers, etc.
 *
 * @class SpriteBatch
 * @extends [PIXI.SpriteBatch](http://www.goodboydigital.com/pixijs/docs/classes/SpriteBatch.html)
 * @uses EventEmitter
 * @constructor
 */
var SpriteBatch = function(settings) {
    PIXI.SpriteBatch.call(this);
    EventEmitter.call(this);

    //mixin user's settings
    utils.setValues(this, settings);

    //Add these properties in so that all objects can see them in the docs
    //these properties are inherited from PIXI.SpriteBatch
    //most of these blocks are copied straight from PIXI source

    /**
     * [read-only] The of children of this object.
     * @property children {Array}
     * @readOnly
     */

    /**
     * The coordinate of the object relative to the local coordinates of the parent.
     *
     * @property position
     * @type Point
     */

    /**
     * The scale factor of the object.
     *
     * @property scale
     * @type Point
     */

    /**
     * The rotation of the object in radians.
     *
     * @property rotation
     * @type Number
     */

    /**
     * The opacity of the object.
     *
     * @property alpha
     * @type Number
     */

    /**
     * The visibility of the object.
     *
     * @property visible
     * @type Boolean
     */

    /**
     * [read-only] The display object that contains this display object.
     *
     * @property parent
     * @type DisplayObject
     * @readOnly
     */

    /**
     * [read-only] The stage the display object is connected to, or undefined if it is not connected to the stage.
     *
     * @property stage
     * @type Stage
     * @readOnly
     */

    /**
     * This is the defined area that will pick up mouse / touch events. It is null by default.
     * Setting it is a neat way of optimising the hitTest function that the interactionManager
     * will use (as it will not need to hit test all the children)
     *
     * @property hitArea
     * @type Rectangle|Polygon|Circle|Ellipse
     */

    /**
     * Wether or not the object will handle mouse events
     *
     * @property interactive
     * @type Boolean
     * @default false
     */
};

inherit(SpriteBatch, PIXI.SpriteBatch, {
    /**
     * Sets the container to visible = true
     *
     * @method show
     * @return {SpriteBatch} Returns itself.
     * @chainable
     */
    show: function() {
        this.visible = true;
        return this;
    },
    /**
     * Sets the container to visible = false
     *
     * @method hide
     * @return {SpriteBatch} Returns itself.
     * @chainable
     */
    hide: function() {
        this.visible = false;
        return this;
    },
    /**
     * Adds a child to the container and returns the child
     *
     * @method addChild
     * @param child {Container|Sprite} Any container or sprite
     * @return {Container|Sprite} The child that was added
     */
    addChild: function(child) {
        PIXI.SpriteBatch.prototype.addChild.apply(this, arguments);

        return child;
    },

    /**
     * Adds a child to the object at a specified index. If the index is out of bounds an error will be thrown
     *
     * @method addChildAt
     * @param child {Container|Sprite} Any container or sprite
     * @param index {Number}
     * @return {Container|Sprite} The child that was added
     */
    addChildAt: function(child) {
        PIXI.SpriteBatch.prototype.addChildAt.apply(this, arguments);

        return child;
    },

    /**
     * Removes a child from the object.
     *
     * @method removeChild
     * @param child {Container|Sprite} Any container or sprite
     * @return {Container|Sprite} The child that was added
     */
    removeChild: function(child) {
        PIXI.SpriteBatch.prototype.removeChild.apply(this, arguments);

        return child;
    },

    /**
     * Removes all children from the object.
     *
     * @method removeAllChildren
     * @return {Container} Returns itself.
     * @chainable
     */
    removeAllChildren: function() {
        while(this.children.length) {
            this.removeChild(this.children[0]);
        }

        return this;
    },

    /**
     * Brings a child to the top of the Z pile.
     *
     * @method bringChildToTop
     * @param child {Container|Sprite} Any container or sprite
     * @return {Container|Sprite} The child that was added
     */
    bringChildToTop: function(child) {
        if(child.parent === this) {
            this.addChild(this.removeChild(child));
        }

        return child;
    },

    /**
     * Destroys this object.
     *
     * @method destroy
     */
    destroy: function() {
        this.disablePhysics();
        this.destroyAllChildren();

        if(this.parent)
            this.parent.removeChild(this);
    },

    /**
     * Destroys all the children of the object.
     *
     * @method destroyAllChildren
     * @return {Container} Returns itself.
     * @chainable
     */
    destroyAllChildren: function() {
        while(this.children.length) {
            if(this.children[0].destroy) {
                this.children[0].destroy();
            } else {
                this.removeChild(this.children[0]);
            }
        }

        return this;
    }
});

module.exports = SpriteBatch;

},{"../utils/EventEmitter":64,"../utils/inherit":69,"../utils/utils":71,"pixi.js":6}],21:[function(_dereq_,module,exports){
/**
 * A texture stores the information that represents an image or part of an image. It cannot be added
 * to the display list directly. It is used to desribe how a Sprite looks. If no frame is provided
 * then the whole image is used.
 * *This is directly exposing [PIXI.Texture](http://www.goodboydigital.com/pixijs/docs/classes/Texture.html)*
 * though it does add some extra methods.
 *
 * @class Texture
 * @uses EventEmitter
 * @constructor
 * @param baseTexture {BaseTexture} The base texture source to create the texture from
 * @param frame {Rectangle} The rectangle frame of the texture to show
 */
var Texture = _dereq_('pixi.js').Texture,
    Rectangle = _dereq_('../geom/Rectangle'),
    utils = _dereq_('../utils/utils'),
    PIXI = _dereq_('pixi.js');

//These create arrays of textures based on texture atlas data

// JSON
Texture.fromJSON = function(key, json, baseTexture) {
    if(!json.frames) {
        utils.warn('Invalid Texture Atlas JSON for fromJSON, missing "frames" array, full json:', json);
        return;
    }

    var frames = json.frames,
        textures = {},
        subkey;

    if(frames.length) {
        for(var i = 0, il = frames.length; i < il; ++i) {
            subkey = key + '_' + frames[i].filename;
            textures[frames[i].filename] = Texture._createFrame(subkey, frames[i], baseTexture);
        }
    } else {
        for(var k in frames) {
            subkey = key + '_' + k;
            textures[k] = Texture._createFrame(subkey, frames[k], baseTexture);
        }
    }

    return textures;
};

Texture._createFrame = function(key, data, baseTexture) {
    var rect = data.frame;

    if(rect) {
        var tx = PIXI.TextureCache[key] = new Texture(baseTexture, {
            x: rect.x,
            y: rect.y,
            width: rect.w,
            height: rect.h
        });

        if(tx.trimmed = data.trimmed) {
            tx.trim = new Rectangle(data.spriteSourceSize.x, data.spriteSourceSize.y, data.sourceSize.w, data.sourceSize.h);
        }

        return tx;
    }
};

// XML
Texture.fromXML = function(key, xml, baseTexture) {
    if(!xml.getElementsByTagName('TextureAtlas')) {
        utils.warn('Invalid Texture Atlas XML given, missing <TextureAtlas> tag, full xml:', xml);
        return;
    }

    var frames = xml.getElementsByTagName('SubTexture') || xml.getElementsByTagName('sprite'),
        textures = {};

    for(var i = 0; i < frames.length; i++) {
        var frame = frames[i],
            attrs = frame.attributes,
            name = attrs.getNamedItem('name') || attrs.getNamedItem('n'),
            //sprite
            x = attrs.getNamedItem('x'),
            y = attrs.getNamedItem('y'),
            width = attrs.getNamedItem('width') || attrs.getNamedItem('w'),
            height = attrs.getNamedItem('height') || attrs.getNamedItem('h'),
            //trim
            ox = attrs.getNamedItem('frameX') || attrs.getNamedItem('oX'),
            oy = attrs.getNamedItem('frameY') || attrs.getNamedItem('oY'),
            owidth = attrs.getNamedItem('frameWidth') || attrs.getNamedItem('oW'),
            oheight = attrs.getNamedItem('frameHeight') || attrs.getNamedItem('oH'),
            //rotated (generic xml export)
            rotated = !!attrs.getNamedItem('r');

        var tx = textures[name] = PIXI.TextureCache[key + '_' + name] = new Texture(baseTexture, {
            x: parseInt(x.nodeValue, 10),
            y: parseInt(y.nodeValue, 10),
            width: parseInt(width.nodeValue, 10),
            height: parseInt(height.nodeValue, 10)
        });

        tx.trimmed = ox && oy;
        tx.rotated = rotated;

        if(tx.trimmed) {
            tx.sourceSize = {
                w: parseInt(owidth.nodeValue, 10),
                h: parseInt(oheight.nodeValue, 10)
            };

            tx.realSize = {
                x: Math.abs(parseInt(ox.nodeValue, 10)),
                y: Math.abs(parseInt(oy.nodeValue, 10)),
                w: parseInt(owidth.nodeValue, 10),
                h: parseInt(oheight.nodeValue, 10)
            };
        }
    }

    return textures;
};

// create textures from sprite sheet
//obj.key
//obj.texture
//obj.image
//obj.frameWidth
//obj.frameHeight
//obj.numFrames
Texture.fromSpritesheet = function(obj) {
    return obj;
};

module.exports = Texture;

},{"../geom/Rectangle":37,"../utils/utils":71,"pixi.js":6}],22:[function(_dereq_,module,exports){
/**
 * A tiling sprite is a fast way of rendering a tiling image
 * *This is directly exposing [PIXI.TilingSprite](http://www.goodboydigital.com/pixijs/docs/classes/TilingSprite.html)*
 *
 * @class TilingSprite
 * @extends [PIXI.DisplayObjectContainer](http://www.goodboydigital.com/pixijs/docs/classes/DisplayObjectContainer.html)
 * @constructor
 * @param texture {Texture} the texture of the tiling sprite
 * @param width {Number}  the width of the tiling sprite
 * @param height {Number} the height of the tiling sprite
 */
var TilingSprite = _dereq_('pixi.js').TilingSprite;

module.exports = TilingSprite;

},{"pixi.js":6}],23:[function(_dereq_,module,exports){
var Effect = _dereq_('./Effect'),
    inherit = _dereq_('../../utils/inherit');

/**
 * Close camera effect. This effect creates a mask on the world that will animated to cover
 * the screen working from the outside-in. It is like a camera shutter "closing" around the target
 *
 * @class fx.camera.Close
 * @extends fx.camera.Effect
 * @constructor
 */
var Close = function() {
    Effect.call(this);
};

inherit(Close, Effect, {
    /**
     * Starts running the effect
     *
     * @method start
     * @param [shape='circle'] {String} The shape to close with, can be either 'ellipse', 'circle', or 'rectangle'
     * @param [duration=1000] {Number} Number of milliseconds for the animation to complete
     * @param [position] {Vector} The position for the animation to close in on, defaults to camera center
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Close} Returns itself.
     * @chainable
     */
    start: function(shape, duration, pos, cb) {
        if(typeof pos ==='function') {
            cb = pos;
            pos = null;
        }

        if(typeof duration === 'function') {
            cb = duration;
            pos = null;
            duration = null;
        }

        if(typeof shape === 'function') {
            cb = shape;
            pos = null;
            duration = null;
            shape = null;
        }

        Effect.prototype.start.call(this, cb);

        this.shape = shape || 'circle';
        this.duration = duration && duration > 0 ? duration : 1000;

        this.cx = pos ? pos.x : this.parent.size.x / 2;
        this.cy = pos ? pos.y : this.parent.size.y / 2;
        this.w = this.mx = this.parent.size.x;
        this.h = this.my = this.parent.size.y;
        this.radius = this.maxRadius = Math.max(this.w / 2, this.h / 2);

        this.gfx.visible = true;
        this.gfx.position.x = this.cx;
        this.gfx.position.y = this.cy;

        this.parent.state.mask = this.gfx;

        if(shape === 'ellipse') {
            this.gfx.scale.y = 0.5;
        }
        else {
            this.gfx.scale.y = 1;
        }

        return this;
    },
    /**
     * Stops running the effect, and removes it from display
     *
     * @method stop
     * @return {fx.camera.Close} Returns itself.
     * @chainable
     */
    stop: function() {
        Effect.prototype.stop.call(this);

        this.radius = this.sx = this.sy = 0;
        this.gfx.visible = false;

        if(this.parent.state.mask === this.gfx)
            this.parent.state.mask = null;

        return this;
    },
    /**
     * Called internally by the camera each frame to update the effect
     *
     * @method update
     * @return {fx.camera.Close} Returns itself.
     * @chainable
     * @private
     */
    update: function(dt) {
        if(this.done) return;

        var part = (dt * 1000) / this.duration;

        this.gfx.clear();
        this.gfx.beginFill(0xff00ff);

        switch(this.shape) {
            case 'ellipse':
            case 'circle':
                this.radius -= (part * this.maxRadius);

                if(this.radius <= 0) {
                    this._complete();
                } else {
                    this.gfx.drawCircle(0, 0, this.radius);
                }
                break;

            case 'rect':
            case 'rectangle':
                this.w -= (part * this.mx);
                this.h -= (part * this.my);

                if(this.w <= 0) {
                    this._complete();
                } else {
                    this.gfx.drawRect(-(this.w / 2), -(this.h / 2), this.w, this.h);
                }
                break;
        }
        this.gfx.endFill();

        return this;
    }
});

module.exports = Close;

},{"../../utils/inherit":69,"./Effect":24}],24:[function(_dereq_,module,exports){
var Container = _dereq_('../../display/Container'),
    inherit = _dereq_('../../utils/inherit'),
    Graphics = _dereq_('../../display/Graphics');

/**
 * Base camera effect class.
 *
 * @class fx.camera.Effect
 * @extends Container
 * @constructor
 */
var Effect = function() {
    Container.call(this);

    /**
     * A graphics instance that can be used by effects to draw
     *
     * @property gfx
     * @type Graphics
     */
    this.gfx = this.addChild(new Graphics());
    this.gfx.visible = false;

    /**
     * Whether or not the effect has completed, and is no longer runnning.
     *
     * @property done
     * @type Boolean
     * @default false
     * @readOnly
     */
    this.done = true;
};

inherit(Effect, Container, {
    /**
     * Starts running the effect
     *
     * @method start
     * @param [callback] {Function} Called when the animation completes
     * @return {Effect} Returns itself.
     * @chainable
     */
    start: function(cb) {
        this.done = false;
        this.cb = cb;
        return this;
    },
    /**
     * Stops running the effect
     *
     * @method stop
     * @return {Effect} Returns itself.
     * @chainable
     */
    stop: function() {
        this.done = true;
        return this;
    },
    /**
     * Called internally by the camera each frame to update the effect
     *
     * @method update
     * @return {Effect} Returns itself.
     * @chainable
     * @private
     */
    update: function() {
        return this;
    },
    /**
     * Called when the effect finishes to call the registered callback (if there is one).
     * If the callback explicitly returns `false` then `.stop()` will not be called. `done`
     * will still be set to `true`but the effect will not be removed from the display.
     *
     * This is useful if you want to run an animation and keep the final state active until
     * you manually remove the item with `.stop()`. For example: fading to black then running
     * some async process, then removing the black manually.
     *
     * @method _complete
     * @private
     */
    _complete: function() {
        this.done = true;

        if(typeof this.cb === 'function') {
            var ret = this.cb();

            if(ret !== false)
                this.stop();
        } else {
            this.stop();
        }
    }
});

module.exports = Effect;

},{"../../display/Container":16,"../../display/Graphics":17,"../../utils/inherit":69}],25:[function(_dereq_,module,exports){
var Effect = _dereq_('./Effect'),
    inherit = _dereq_('../../utils/inherit');

/**
 * Fade the screen into a color. This will fade into a color that will
 * eventually cover the screen.
 *
 * @class fx.camera.Fade
 * @extends fx.camera.Effect
 * @constructor
 */
var Fade = function() {
    Effect.call(this);
};

inherit(Fade, Effect, {
    /**
     * Starts running the effect
     *
     * @method start
     * @param [color=0xFFFFFF] {Number} The color to fade into
     * @param [duration=1000] {Number} The time it should take (in milliseconds) to fade in
     * @param [alpha=1] {Number} The opacity to fade into (final opacity)
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Fade} Returns itself.
     * @chainable
     */
    start: function(color, duration, alpha, cb) {
        if(typeof alpha === 'function') {
            cb = duration;
            alpha = null;
        }

        if(typeof duration === 'function') {
            cb = duration;
            alpha = null;
            duration = null;
        }

        if(typeof color === 'function') {
            cb = color;
            alpha = null;
            duration = null;
            color = null;
        }

        Effect.prototype.start.call(this, cb);

        color = typeof color === 'number' ? color : 0xFFFFFF;
        this.goal = alpha || 1;
        this.duration = duration && duration > 0 ? duration : 1000;

        this.gfx.visible = true;
        this.gfx.alpha = 0;
        this.gfx.clear();
        this.gfx.beginFill(color);
        this.gfx.drawRect(0, 0, this.parent.size.x, this.parent.size.y);

        return this;
    },
    /**
     * Stops running the effect, and removes it from display
     *
     * @method stop
     * @return {fx.camera.Fade} Returns itself.
     * @chainable
     */
    stop: function() {
        Effect.prototype.stop.call(this);

        this.gfx.alpha = 0;
        this.gfx.visible = false;

        return this;
    },
    /**
     * Called internally by the camera each frame to update the effect
     *
     * @method update
     * @return {fx.camera.Fade} Returns itself.
     * @chainable
     * @private
     */
    update: function(dt) {
        if(this.done) return;

        if(this.gfx.alpha < this.goal) {
            this.gfx.alpha += (dt * 1000) / this.duration;

            if(this.gfx.alpha >= this.goal) {
                this._complete();
            }
        }

        return this;
    }
});

module.exports = Fade;

},{"../../utils/inherit":69,"./Effect":24}],26:[function(_dereq_,module,exports){
var Effect = _dereq_('./Effect'),
    inherit = _dereq_('../../utils/inherit');

/**
 * Flash the screen with a color. This will cover the screen in a
 * color then fade it out.
 *
 * @class fx.camera.Flash
 * @extends fx.camera.Effect
 * @constructor
 */
var Flash = function() {
    Effect.call(this);
};

inherit(Flash, Effect, {
    /**
     * Starts running the effect
     *
     * @method start
     * @param [color=0xFFFFFF] {Number} The color to flash the screen with
     * @param [duration=1000] {Number} The time it should take (in milliseconds) to fade out
     * @param [alpha=1] {Number} The opacity of the initial flash of color (start opacity)
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Flash} Returns itself.
     * @chainable
     */
    start: function(color, duration, alpha, cb) {
        if(typeof alpha === 'function') {
            cb = duration;
            alpha = null;
        }

        if(typeof duration === 'function') {
            cb = duration;
            alpha = null;
            duration = null;
        }

        if(typeof color === 'function') {
            cb = color;
            alpha = null;
            duration = null;
            color = null;
        }

        Effect.prototype.start.call(this, cb);

        alpha = alpha || 1;
        color = typeof color === 'number' ? color : 0xFFFFFF;
        this.duration = duration && duration > 0 ? duration : 1000;

        this.gfx.visible = true;
        this.gfx.alpha = alpha;
        this.gfx.clear();
        this.gfx.beginFill(color);
        this.gfx.drawRect(0, 0, this.parent.size.x, this.parent.size.y);

        return this;
    },
    /**
     * Stops running the effect, and removes it from display
     *
     * @method stop
     * @return {fx.camera.Flash} Returns itself.
     * @chainable
     */
    stop: function() {
        Effect.prototype.stop.call(this);

        this.gfx.alpha = 0;
        this.gfx.visible = false;

        return this;
    },
    /**
     * Called internally by the camera each frame to update the effect
     *
     * @method update
     * @return {fx.camera.Flash} Returns itself.
     * @chainable
     * @private
     */
    update: function(dt) {
        if(this.done) return;

        if(this.gfx.alpha > 0) {
            this.gfx.alpha -= (dt * 1000) / this.duration;

            if(this.gfx.alpha <= 0) {
                this._complete();
            }
        }

        return this;
    }
});

module.exports = Flash;

},{"../../utils/inherit":69,"./Effect":24}],27:[function(_dereq_,module,exports){
var Effect = _dereq_('./Effect'),
    inherit = _dereq_('../../utils/inherit'),
    C = _dereq_('../../constants');

/**
 * Adds arcade-style scanlines to the camera viewport.
 *
 * @class fx.camera.Scanlines
 * @extends fx.camera.Effect
 * @constructor
 */
var Scanlines = function() {
    Effect.call(this);
};

inherit(Scanlines, Effect, {
    /**
     * Starts running the effect
     *
     * @method start
     * @param [color=0x000000] {Number} The color for the scanlines to be
     * @param [axis=gf.AXIS.HORIZONTAL] {gf.AXIS} The axis to draw the lines on
     * @param [spacing=4] {Number} Number of pixels between each line
     * @param [thickness=1] {Number} Number of pixels thick each line is
     * @param [alpha=0.3] {Number} The opacity of the lines
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Scanlines} Returns itself.
     * @chainable
     */
    start: function(color, axis, spacing, thickness, alpha, cb) {
        if(typeof alpha ==='function') {
            cb = alpha;
            alpha = null;
        }

        if(typeof thickness === 'function') {
            cb = thickness;
            alpha = null;
            thickness = null;
        }

        if(typeof spacing === 'function') {
            cb = spacing;
            alpha = null;
            thickness = null;
            spacing = null;
        }

        if(typeof axis === 'function') {
            cb = spacing;
            alpha = null;
            thickness = null;
            spacing = null;
            axis = null;
        }

        if(typeof color === 'function') {
            cb = spacing;
            alpha = null;
            thickness = null;
            spacing = null;
            axis = null;
            color = null;
        }

        Effect.prototype.start.call(this, cb);

        color = color || 0x000000;
        axis = axis || C.AXIS.HORIZONTAL;
        spacing = spacing || 4;
        thickness = thickness || 1;
        alpha = alpha || 0.3;

        var sx = this.parent.size.x,
            sy = this.parent.size.y;

        this.gfx.clear();
        this.gfx.visible = true;
        this.gfx.beginFill(color, alpha);

        //draw the lines
        if(axis & C.AXIS.VERTICAL) {
            for(var x = 0; x < sx; x += spacing) {
                this.gfx.drawRect(x, 0, thickness, sy);
            }
        }

        if(axis & C.AXIS.HORIZONTAL) {
            for(var y = 0; y < sy; y += spacing) {
                this.gfx.drawRect(0, y, sx, thickness);
            }
        }
        this.gfx.endFill();

        return this;
    },
    /**
     * Stops running the effect, and removes it from display
     *
     * @method stop
     * @return {fx.camera.Scanlines} Returns itself.
     * @chainable
     */
    stop: function() {
        Effect.prototype.stop.call(this);

        this.gfx.visible = false;

        return this;
    }
});

module.exports = Scanlines;

},{"../../constants":11,"../../utils/inherit":69,"./Effect":24}],28:[function(_dereq_,module,exports){
var Effect = _dereq_('./Effect'),
    Vector = _dereq_('../../math/Vector'),
    inherit = _dereq_('../../utils/inherit'),
    math = _dereq_('../../math/math'),
    C = _dereq_('../../constants');

/**
 * Shakes the camera around a bit.
 *
 * @class fx.camera.Shake
 * @extends fx.camera.Effect
 * @constructor
 */
var Shake = function() {
    Effect.call(this);
    this.offset = new Vector();
};

inherit(Shake, Effect, {
    /**
     * Starts running the effect
     *
     * @method start
     * @param [intensity=0.01] {Number} The intensity of the shaking
     * @param [duration=1000] {Number} The amount of time the screen shakes for (in milliseconds)
     * @param [direction=gf.AXIS.BOTH] {gf.AXIS} The axis to shake on
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Shake} Returns itself.
     * @chainable
     */
    start: function(intensity, duration, direction, cb) {
        if(typeof direction === 'function') {
            cb = direction;
            direction = null;
        }

        if(typeof duration === 'function') {
            cb = duration;
            direction = null;
            duration = null;
        }

        if(typeof intensity === 'function') {
            cb = intensity;
            direction = null;
            duration = null;
            intensity = null;
        }

        Effect.prototype.start.call(this, cb);

        this.intensity = intensity || 0.01;
        this.duration = duration || 1000;
        this.direction = direction || C.AXIS.BOTH;
        this.offset.x = this.offset.y = 0;

        return this;
    },
    /**
     * Stops running the effect, and removes it from display
     *
     * @method stop
     * @return {fx.camera.Shake} Returns itself.
     * @chainable
     */
    stop: function() {
        Effect.prototype.stop.call(this);

        this.duration = this.offset.x = this.offset.y = 0;

        return this;
    },
    /**
     * Called internally by the camera each frame to update the effect
     *
     * @method update
     * @return {fx.camera.Shake} Returns itself.
     * @chainable
     * @private
     */
    update: function(dt) {
        if(this.done) return;

        this.duration -= (dt * 1000);

        //pan back to the original position
        this.offset.x = -this.offset.x;
        this.offset.y = -this.offset.y;
        this.parent.pan(this.offset.x, this.offset.y);

        //check if we are complete
        if(this.duration <= 0) {
            this._complete();
        }
        //otherwise do the shake
        else {
            //pan to a random offset
            if(this.direction & C.AXIS.HORIZONTAL)
                this.offset.x = math.round(Math.random() * this.intensity * this.parent.size.x * 2 - this.intensity * this.parent.size.x);

            if (this.direction & C.AXIS.VERTICAL)
                this.offset.y = math.round(Math.random() * this.intensity * this.parent.size.y * 2 - this.intensity * this.parent.size.y);

            this.parent.pan(this.offset.x, this.offset.y);
        }
    }
});

module.exports = Shake;

},{"../../constants":11,"../../math/Vector":47,"../../math/math":48,"../../utils/inherit":69,"./Effect":24}],29:[function(_dereq_,module,exports){
var PIXI = _dereq_('pixi.js');

/**
 *
 * @class Filter
 * @extends PIXI.AbstractFilter
 * @constructor
 * @param fragmentSrc {Array<String>} The source of the fragment shader
 * @param uniforms {Object} The uniforms for this filter
 */
var Filter = PIXI.AbstractFilter;

/**
 * An array of passes - some filters contain a few steps this array simply stores the steps in a liniear fashion.
 * For example the blur filter has two passes blurX and blurY.
 *
 * @property passes
 * @type An array of filter objects
 * @private
 */

/**
 * The raw PIXI.Shader instances that were created by a WebGLRenderer.
 *
 * @property shaders
 * @type An array of shader objects
 * @private
 */

/**
 * The uniforms used by this filter (and associated fragment shader)
 *
 * @property uniforms
 * @type Object
 * @private
 */

/**
 * The fragment shader source in a line-by-line array
 *
 * @property fragmentSrc
 * @type Array<String>
 * @private
 */

module.exports = Filter;

},{"pixi.js":6}],30:[function(_dereq_,module,exports){
var StateManager = _dereq_('./StateManager'),
    EventEmitter = _dereq_('../utils/EventEmitter'),
    Cache = _dereq_('../utils/Cache'),
    Clock = _dereq_('../utils/Clock'),
    SpritePool = _dereq_('../utils/SpritePool'),
    Loader = _dereq_('../loader/Loader'),
    InputManager = _dereq_('../input/InputManager'),
    AudioManager = _dereq_('../audio/AudioManager'),
    Vector = _dereq_('../math/Vector'),
    utils = _dereq_('../utils/utils'),
    support = _dereq_('../utils/support'),
    inherit = _dereq_('../utils/inherit'),
    PIXI = _dereq_('pixi.js'),
    C = _dereq_('../constants');

/**
 * Main game object, controls the entire instance of the game
 *
 * @class Game
 * @extends Object
 * @uses EventEmitter
 * @constructor
 * @param container {DOMElement|String} The container for the new canvas we will create for the game, or the ID of one
 * @param settings {Object} All the settings for the game instance
 * @param settings.width {Number} The width of the viewport
 * @param settings.height {Number} The height of the viewport
 * @param [settings.renderer=RENDERER.AUTO] {String} The renderer to use either RENDERER.AUTO, RENDERER.CANVAS, or RENDERER.WEBGL
 * @param [settings.transparent=false] {Boolean} Should the render element have a transparent background
 * @param [settings.background='#FFF'] {Number} The background color of the stage
 * @param [settings.antialias=true] {Boolean} Anti-alias graphics (in WebGL this helps with edges, in Canvas2D it retains pixel-art quality)
 * @param [settings.canvas] {DOMElement} The canvas to render into, if not specified one is created
 */
var Game = function(container, settings) {
    EventEmitter.call(this);

    //setup settings defaults
    settings = settings || {};
    settings.width = settings.width || 800;
    settings.height = settings.height || 600;
    settings.renderer = settings.renderer || C.RENDERER.AUTO;
    settings.transparent = settings.transparent || false;
    settings.background = settings.background || '#FFF';
    settings.antialias = settings.antialias !== undefined ? settings.antialias : true;
    settings.canvas = settings.canvas || null; //passing null to renderer, lets the renderer make one

    /**
     * The domElement that we are putting our rendering canvas into (the container)
     *
     * @property container
     * @type DOMELement
     * @readOnly
     */
    this.container = typeof container === 'string' ? document.getElementById(container) : container;

    if(!this.container)
        this.container = document.body;

    /**
     * The width of the render viewport
     *
     * @property width
     * @type Number
     * @default 800
     */
    this.width = settings.width;

    /**
     * The height of the render viewport
     *
     * @property height
     * @type Number
     * @default 600
     */
    this.height = settings.height;

    /**
     * The method used to render values to the screen (either webgl, or canvas)
     *
     * @property renderMethod
     * @type String
     * @default RENDERER.AUTO
     */
    this.renderMethod = settings.renderer;

    /**
     * Whether the canvas has a transparent background or not
     *
     * @property transparent
     * @type Boolean
     * @default false
     */
    this.transparent = settings.transparent;

    /**
     * The background of the stage
     *
     * @property background
     * @type Boolean
     * @default false
     */
    this.background = settings.background;

    /**
     * Anti-alias graphics (in WebGL this helps with edges, in Canvas2D it retains pixel-art quality)
     *
     * @property antialias
     * @type Boolean
     * @default true
     */
    this.antialias = settings.antialias;

    /**
     * The canvas to render into
     *
     * @property canvas
     * @type HTMLCanvasElement
     */
    this.canvas = settings.canvas;

    /**
     * Raw rendering engine, the underlying PIXI renderer that draws for us
     *
     * @property renderer
     * @type PIXI.WebGLRenderer|PIXI.CanvasRenderer
     * @readOnly
     */
    this.renderer = this._createRenderer();

    /**
     * Raw PIXI.stage instance, the root of all things in the scene graph
     *
     * @property stage
     * @type PIXI.Stage
     * @readOnly
     */
    this.stage = new PIXI.Stage(this.background);

    /**
     * Clock instance for internal timing
     *
     * @property clock
     * @type Clock
     * @readOnly
     */
    this.clock = new Clock();

    /**
     * The audio manager for this game instance, used to play and control
     * all the audio in a game.
     *
     * @property audio
     * @type AudioManager
     * @readOnly
     */
    this.audio = new AudioManager(this);

    /**
     * The loader for this game instance, used to preload assets into the cache
     *
     * @property loader
     * @type Loader
     * @readOnly
     */
    this.load = new Loader(this);

    /**
     * Cache instance for storing/retrieving assets
     *
     * @property cache
     * @type Cache
     * @readOnly
     */
    this.cache = new Cache(this);

    /**
     * The input instance for this game
     *
     * @property input
     * @type InputManager
     * @readOnly
     */
    this.input = new InputManager(this);

    /**
     * The sprite pool to use to create registered entities
     *
     * @property spritepool
     * @type SpritePool
     * @readOnly
     */
    this.spritepool = new SpritePool(this);

    /**
     * The state manager, to switch between game states
     *
     * @property state
     * @type StateManager
     * @readOnly
     */
    this.state = new StateManager(this);

    /**
     * The offset for the viewport in the document
     *
     * @property offset
     * @type Vector
     * @readOnly
     */
    this.offset = new Vector();

    /**
     * Holds timing data for the previous loop
     *
     * @property timings
     * @type Object
     * @readOnly
     */
    this.timings = {};

    this.maxDelta = 5;

    //pixi does some prevent default on mousedown, so we need to
    //make sure mousedown will focus the canvas or keyboard events break
    var view = this.canvas;

    if(!view.getAttribute('tabindex'))
        view.setAttribute('tabindex','1');

    view.focus();
    view.addEventListener('click', function() {
        view.focus();
    }, false);

    /**
     * Fired each frame after everything has updated, but just before rendering
     *
     * @event tick
     * @param dt {Number} The number of seconds passed since the last tick call (delta time)
     */
};

inherit(Game, Object, {
    /**
     * Creates the underlying renderer based on browser support. It will also set's `game.renderMethod` for a user
     * to be able to check.
     *
     * @method _createRenderer
     * @return {PIXI.WebGLRenderer|PIXI.CanvasRenderer} The renderer to use
     * @private
     */
    _createRenderer: function() {
        var method = this.renderMethod,
            render = null;

        //no support
        if(!support.webgl && !support.canvas) {
            throw new Error('Neither WebGL nor Canvas is supported by this browser!');
        }
        else if((method === C.RENDERER.WEBGL || method === C.RENDERER.AUTO) && support.webgl) {
            method = C.RENDERER.WEBGL;
            render = new PIXI.WebGLRenderer(this.width, this.height, this.canvas, this.transparent, this.antialias);
        }
        else if((method === C.RENDERER.CANVAS || method === C.RENDERER.AUTO) && support.canvas) {
            method = C.RENDERER.CANVAS;
            render = new PIXI.CanvasRenderer(this.width, this.height, this.canvas, this.transparent);
        }
        else {
            throw new Error('Your render method ("' + method + '") is not supported by this browser!');
        }

        //append the renderer view only if the user didn't pass their own
        if(!this.canvas) {
            this.container.appendChild(render.view);
            this.canvas = render.view;
        }

        this.offset = utils.getOffset(this.canvas);
        this.renderMethod = method;

        return render;
    },
    /**
     * Allows you to resize the game area.
     *
     * @method resize
     * @param width {Number} Width to resize to
     * @param height {Number} Height to resize to
     * @return {Game} Returns itself.
     * @chainable
     */
    resize: function(w, h) {
        this.renderer.resize(w, h);
        this.width = w;
        this.height = h;

        for(var i = 0, il = this.stage.children.length; i < il; ++i) {
            var o = this.stage.children[i];

            if(o.resize)
                o.resize(w, h);
        }

        return this;
    },
    /**
     * Requests that the browser go into fullscreen mode.
     *
     * @method requestFullscreen
     * @return {Game} Returns itself.
     * @chainable
     */
    requestFullscreen: function() {
        var elem = this.renderer.view;

        if(elem.requestFullscreen) {
            elem.requestFullscreen();
        } else if(elem.mozRequestFullScreen) {
            elem.mozRequestFullScreen();
        } else if(elem.webkitRequestFullscreen) {
            elem.webkitRequestFullscreen();
        }

        return this;
    },
    /**
     * Begins the render loop.
     *
     * @method render
     * @return {Game} Returns itself.
     * @chainable
     */
    render: function() {
        this.clock.start();
        this._tick();

        return this;
    },
    /**
     * The looping render tick.
     *
     * @method _tick
     * @private
     */
    _tick: function() {
        //start render loop
        window.requestAnimFrame(this._tick.bind(this));

        var dt = this.clock.getDelta();

        if(dt > this.maxDelta) return;

        this.timings.lastTickStart = this.timings.tickStart;
        this.timings.tickStart = this.clock.now();

        this.timings.lastDelta = dt;

        //gather input from user
        this.timings.inputStart = this.clock.now();
        this.input.update(dt);
        this.timings.inputEnd = this.clock.now();

        //TODO: plugins
        //this.timings.pluginsStart = this.clock.now();
        //this.plugins.update(dt);
        //this.timings.pluginsEnd = this.clock.now();

        //update this game state
        this.timings.stateStart = this.clock.now();
        this.state.active.update(dt);
        this.timings.stateEnd = this.clock.now();

        this.timings.userFuncsStart = this.clock.now();
        this.emit('tick', dt);
        this.timings.userFuncsEnd = this.clock.now();

        //render scene
        this.timings.renderStart = this.clock.now();
        this.renderer.render(this.stage);
        this.timings.renderEnd = this.clock.now();

        this.timings.tickEnd = this.clock.now();
    }
});

/**
 * Alias for the active State's physics object. Instead of using
 * `game.state.active.physics`, you can use `game.physics`
 *
 * @property physics
 * @type Physics
 * @readOnly
 */
Object.defineProperty(Game.prototype, 'physics', {
    get: function() {
        return this.state.active.physics;
    }
});

/**
 * Alias for the active State's camera object. Instead of using
 * `game.state.active.camera`, you can use `game.camera`
 *
 * @property camera
 * @type Camera
 * @readOnly
 */
Object.defineProperty(Game.prototype, 'camera', {
    get: function() {
        return this.state.active.camera;
    }
});

/**
 * Alias for the active State's world object. Instead of using
 * `game.state.active.world`, you can use `game.world`
 *
 * @property world
 * @type World
 * @readOnly
 */
Object.defineProperty(Game.prototype, 'world', {
    get: function() {
        return this.state.active.world;
    }
});

module.exports = Game;

},{"../audio/AudioManager":8,"../constants":11,"../input/InputManager":40,"../loader/Loader":46,"../math/Vector":47,"../utils/Cache":61,"../utils/Clock":62,"../utils/EventEmitter":64,"../utils/SpritePool":68,"../utils/inherit":69,"../utils/support":70,"../utils/utils":71,"./StateManager":32,"pixi.js":6}],31:[function(_dereq_,module,exports){
var AudioManager = _dereq_('../audio/AudioManager'),
    Container = _dereq_('../display/Container'),
    World = _dereq_('./World'),
    Camera = _dereq_('../camera/Camera'),
    PhysicsSystem = _dereq_('../physics/PhysicsSystem'),
    math = _dereq_('../math/math'),
    inherit = _dereq_('../utils/inherit');

/**
 * States are containers that represent different states of a game
 *
 * @class State
 * @extends Container
 * @constructor
 * @param game {Game} The game instance this state belongs to
 * @param [name] {String} The name of this state
 * @param [physicsOptions] {Object} All the settings for the physics environment
 * @param [physicsOptions.gravity] {Vector} The gravity constant for the physics system (default is 9.87, which is normal Earth gravity)
 * @example
 *      var state = game.state.add('battle');
 *      state.addChild(battlePlayer);
 *      state.addChild(enemy);
 *
 *      state.enable();
 */
var State = function(game, name, physOptions) {
    if(!name)
        name = math.randomString();

    /**
     * The name of this game state
     *
     * @property name
     * @type String
     */
    this.name = name;

    /**
     * The game instance that this state belongs too, will be set
     * when setup() is called with a game instance.
     *
     * @property game
     * @type Game
     */
    this.game = game;

    /**
     * The audio manager for this game state
     *
     * @property audio
     * @type AudioManager
     * @readOnly
     */
    this.audio = new AudioManager(game, game.audio);

    /**
     * The container that holds all non-gui sprites and the tilemap
     *
     * @property world
     * @type Tilemap
     * @readOnly
     */
    this.world = new World(this);

    /**
     * The physics system to simulate the world physics
     *
     * @property physics
     * @type Physics
     * @readOnly
     */
    this.physics = new PhysicsSystem(this, physOptions);

    /**
     * The camera you view the scene through, will be set
     * when setup() is called with a game instance.
     *
     * @property camera
     * @type Camera
     * @readOnly
     */
    this.camera = new Camera(this);

    //call base ctor
    Container.call(this);

    //start disabled
    this.visible = false;

    //add world/camera
    this.addChild(this.world);
    this.addChild(this.camera);

    //ensure the camera is the right size
    this.camera.resize(game.width, game.height);
};

inherit(State, Container, {
    /**
     * Enables (shows) the game state
     *
     * @method enable
     * @return {State} Returns itself.
     * @chainable
     */
    enable: function() {
        this.game.state.enable(this);

        return this;
    },
    /**
     * Called internally by `game.resize`. This ensures that the camera
     * is the correct size, and renders the world with the new viewport size.
     *
     * @method resize
     * @param width {Number} The width of the new viewport
     * @param height {Number} The height of the new viewport
     * @return {State} Returns itself.
     * @chainable
     * @private
     */
    resize: function(w, h) {
        this.camera.resize(w, h);
        this.world.resize(w, h);

        return this;
    },
    /**
     * Called by the game each frame to update the state objects.
     *
     * @method update
     * @param dt {Number} The number of seconds passed since the last update call.
     * @private
     */
    update: function(dt) {
        //update any camera effects
        this.game.timings.cameraStart = this.game.clock.now();
        this.camera.update(dt);
        this.game.timings.cameraEnd = this.game.clock.now();

        //simulate physics and detect/resolve collisions
        this.game.timings.physicsStart = this.game.clock.now();
        this.physics.update(dt);
        this.game.timings.physicsEnd = this.game.clock.now();

        return this;
    }
});

module.exports = State;

},{"../audio/AudioManager":8,"../camera/Camera":10,"../display/Container":16,"../math/math":48,"../physics/PhysicsSystem":52,"../utils/inherit":69,"./World":33}],32:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    State = _dereq_('./State');

/**
 * A state manager is a container for all the states in a game.
 *
 * @class StateManager
 * @extends Object
 * @constructor
 * @param game {Game} The game this manager bleongs to.
 */
var StateManager = function(game) {
    /**
     * The game instance that this manager belongs to.
     *
     * @property game
     * @type Game
     */
    this.game = game;

    /**
     * The states managed by this manager, keyed on the state name
     *
     * @property states
     * @type Object<State>
     */
    this.states = {};

    /**
     * The currently active state
     *
     * @property active
     * @type State
     */
    this.active = null;

    /**
     * The count of states in this manager
     *
     * @property count
     * @type Number
     */
    this.count = 0;

    //create a default state
    this._createDefault();
};

inherit(StateManager, Object, {
    /**
     * Creates the default state
     *
     * @method _createDefault
     * @return {State} The default state
     */
    _createDefault: function() {
        return this.add('__default', true);
    },
    /**
     * Adds a state to the game, creating one if necessary.
     *
     * There are 3 ways to use this function to add a state to the manager. The simplest case
     * is to pass a string for the name, and let the manager create a normal gf.State for you
     * with the name you provided. The second usage is to pass a class that is a decendant of gf.State.
     *
     * For example:
     *
     * ```
     * function MyState(game) {
     *     gf.State.call(game, 'some-name');
     * }
     * gf.inherit(MyState, gf.State);
     *
     * game.state.add(MyState); //adds a new instance of your state
     * ```
     *
     * The final usage is to pass a state that is already created. In this case the manager will
     * add the state to the list based on `state.name` and set the game to be the manager's game
     * instance with `state.game = this.game`;
     *
     * @method add
     * @param state {String|Function|State} The state name, constructor, or state instance to add.
     * @return {State} The state that was added
     */
    add: function(Name, enable) {
        var state;

        //create a state if a string is passed
        if(typeof Name === 'string') {
            state = new State(this.game, Name);
        }
        //create a state of the instance passed
        else if(typeof Name === 'function') {
            state = new Name(this.game);
        }
        //a pre-created state, ensure game is set correctly
        else {
            state = Name;
            state.game = this.game;
        }

        this.states[state.name] = state;
        this.game.stage.addChild(state);

        if(enable)
            this.enable(state);

        this.count++;

        return state;
    },
    /**
     * Removes a state from the game
     *
     * @method remove
     * @param state {String|State} The name of the state to remove, or the state instance itself.
     * @return {StateManager} Returns itself.
     * @chainable
     */
    remove: function(state) {
        if(typeof state === 'string')
            state = this.states[state];

        if(state.parent)
            state.parent.removeChild(state);

        delete this.states[state.name];

        this.count--;

        return this;
    },
    /**
     * Enables a state in the game.
     *
     * @method enable
     * @param state {String|State} The name of the state to enable, or the state instance itself.
     * @return {StateManager} Returns itself.
     * @chainable
     */
    enable: function(state) {
        if(typeof state !== 'string')
            state = state.name;

        if(this.states[state]) {
            if(this.active) {
                this.active.visible = false;
                this.active.physics.pause();
            }

            this.active = this.states[state];
            this.active.visible = true;
            this.active.physics.resume();
        }

        return this;
    },
    /**
     * Destroys the state manager completely
     *
     * @method destroy
     */
    destroy: function() {
        this.game = null;
        this.states = null;
    }
});

module.exports = StateManager;

},{"../utils/inherit":69,"./State":31}],33:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    Container = _dereq_('../display/Container'),
    Rectangle = _dereq_('../geom/Rectangle'),
    ObjectFactory = _dereq_('../utils/ObjectFactory'),
    math = _dereq_('../math/math');

/**
 * The world is the container for all game objects in a game state.
 *
 * @class World
 * @extends Container
 * @constructor
 */
var World = function(state) {
    Container.call(this);

    /**
     * The game instance this world belongs to
     *
     * @property game
     * @type Game
     */
    this.game = state.game;

    /**
     * The game state this world belongs to
     *
     * @property state
     * @type State
     */
    this.state = state;

    /**
     * The bounds of the world
     *
     * @property bounds
     * @type Rectangle
     */
    this.bounds = new Rectangle(0, 0, state.game.width, state.game.height);

    /**
     * An object factory for creating game objects
     *
     * @property add
     * @type ObjectFactory
     */
    this.add = new ObjectFactory(state, this);
};

inherit(World, Container, {
    /**
     * Pans the world around
     *
     * @method pan
     * @param x {Number|Point} The x amount to pan, if a Vector is passed the y param is ignored
     * @param y {Number} The y ammount to pan
     * @return {World} Returns itself.
     * @chainable
     */
    pan: function(x, y) {
        y = math.floor(x.y !== undefined ? x.y : (y || 0));
        x = math.floor(x.x !== undefined ? x.x : (x || 0));

        this.position.x += x * this.scale.x;
        this.position.y += y * this.scale.y;

        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.pan)
                o.pan(x, y);
        }

        return this;
    },
    /**
     * Resizes the children of the world, called by game.resize()
     *
     * @method resize
     * @param width {Number} Width to resize to
     * @param height {Number} Height to resize to
     * @return {World} Returns itself.
     * @chainable
     */
    resize: function(w, h) {
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.render)
                o.render(-this.position.x, -this.position.y, w, h);
        }

        return this;
    }
});

module.exports = World;

},{"../display/Container":16,"../geom/Rectangle":37,"../math/math":48,"../utils/ObjectFactory":65,"../utils/inherit":69}],34:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    Vector = _dereq_('../math/Vector'),
    C = _dereq_('../constants');

/**
 * The Circle object is an area defined by its position, as indicated by its
 * center point (x, y) and by its radius.
 *
 * @class Circle
 * @constructor
 * @param center {Vector} The point of the center of the circle
 * @param radius {Number} The radius of the circle
 */
var Circle = function(x, y, radius, scale) {
    /**
     * The center of the circle
     *
     * @property position
     * @type Vector
     */
    this.position = new Vector();

    /**
     * The unscaled radius of the circle
     *
     * @property _radius
     * @type Number
     * @default 0
     * @private
     */
    this._radius = radius || 0;

    /**
     * The radius of the circle
     *
     * @property radius
     * @type Number
     * @default 0
     */
    this.radius = radius || 0;

    /**
     * The scale of the circle
     *
     * @property scale
     * @type Vector
     * @default new Vector(1, 1)
     */
    this.scale = scale || new Vector(1, 1);

    //set position
    this.x = x || 0;
    this.y = y || 0;

    //internal shape type
    this._shapetype = C.SHAPE.CIRCLE;

    this.recalc();
};

inherit(Circle, Object, {
    /**
     * Creates a clone of this Circle instance
     *
     * @method clone
     * @return {Circle} a copy of the circle
     */
    clone: function() {
        return new Circle(this.x, this.y, this.radius);
    },

    /**
     * Copies the values from another circle to this one
     *
     * @method copy
     * @param circle {Circle} The circle to copy vlaues from
     * @return {Circle} Returns itself.
     * @chainable
     */
    copy: function(circle) {
        this.x = circle.x;
        this.y = circle.y;
        this.radius = circle.radius;

        return this;
    },

    /**
     * Checks if the x, and y coords passed to this function are contained within this circle,
     * or on the edge of the circle
     *
     * @method contains
     * @param x {Number} The X coord of the point to test
     * @param y {Number} The Y coord of the point to test
     * @return {Boolean} if the x/y coords are within this polygon
     */
    contains: function(x, y) {
        if(this.radius <= 0)
            return false;

        var dx = (x - this.x),
            dy = (y - this.y),
            r2 = this.radius * this.radius;

        dx *= dx;
        dy *= dy;

        return (dx + dy <= r2);
    },

    /**
     * Checks if this circle overlaps another
     *
     * @method overlaps
     * @param circle {Circle} The circle to check if this overlaps
     * @return {Boolean} if the circle overlaps
     */
    overlaps: function(circle) {
        var differenceV = this.position.clone().sub(circle.position),
            totalRadius = this.radius + circle.radius,
            totalRadiusSq = totalRadius * totalRadius,
            distanceSq = differenceV.lengthSq();

        //if distanceSq is greater than totalRadiusSq then they do not intersect,
        //so we return the inverse of that value.
        /*jshint -W018*/
        return !(distanceSq > totalRadiusSq);
    },

    /**
     * Checks if this circle's values are equal to anothers
     *
     * @method equals
     * @param circle {Circle} The circle to check against
     * @return {Boolean} True if they are equal
     */
    equals: function(circle) {
        return this.position.equals(circle.position) &&
                this.radius === circle.radius;
    },

    /**
     * Recalculates the scaled radius
     *
     * @method recalc
     * @return {Circle} Returns itself.
     * @chainable
     */
    recalc: function() {
        this.radius = this._radius * this.scale.x;

        return this;
    }
});

/**
 * The center X coord of the circle
 *
 * @property x
 * @type Number
 * @default 0
 */
Object.defineProperty(Circle.prototype, 'x', {
    get: function() {
        return this.position.x;
    },
    set: function(v) {
        this.position.x = v;
    }
});

/**
 * The center Y coord of the circle
 *
 * @property y
 * @type Number
 * @default 0
 */
Object.defineProperty(Circle.prototype, 'y', {
    get: function() {
        return this.position.y;
    },
    set: function(v) {
        this.position.y = v;
    }
});

/**
 * The radius circle
 *
 * @property radius
 * @type Number
 * @default 0
 */
Object.defineProperty(Circle.prototype, 'radius', {
    get: function() {
        return this._radius * this.scale.x;
    },
    set: function(v) {
        this._radius = v;
    }
});

/**
 * The circumference of the circle
 *
 * @property circumference
 * @type Number
 * @readOnly
 */
Object.defineProperty(Circle.prototype, 'circumference', {
    get: function() {
        return 2 * (Math.PI * this.radius);
    }
});

/**
 * The area of the circle
 *
 * @property area
 * @type Number
 * @readOnly
 */
Object.defineProperty(Circle.prototype, 'area', {
    get: function() {
        return Math.PI * this.radius * this.radius;
    }
});

module.exports = Circle;

},{"../constants":11,"../math/Vector":47,"../utils/inherit":69}],35:[function(_dereq_,module,exports){
/**
 * The Ellipse object can be used to specify a hit area for displayobjects
 * see [PIXI.Ellipse](http://www.goodboydigital.com/pixijs/docs/classes/Ellipse.html)
 * for more information.
 *
 * @class Ellipse
 * @constructor
 * @param x {Number} The X coord of the upper-left corner of the framing rectangle of this ellipse
 * @param y {Number} The Y coord of the upper-left corner of the framing rectangle of this ellipse
 * @param width {Number} The overall height of this ellipse
 * @param height {Number} The overall width of this ellipse
 */
var Ellipse = _dereq_('pixi.js').Ellipse;

module.exports = Ellipse;

},{"pixi.js":6}],36:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    Vector = _dereq_('../math/Vector'),
    C = _dereq_('../constants');

/**
 * A *convex* clockwise Polygon.
 *
 * @class Polygon
 * @constructor
 * @param x {Number} The X origin of the polygon, all X coords for all points are relative to this
 * @param y {Number} The Y origin of the polygon, all Y coords for all points are relative to this
 * @param points {Array<Vector>|Array<Number>} This can be an array of Vectors that form the polygon,
 *      a flat array of numbers that will be interpreted as [x,y, x,y, ...]
 * @param scale {Number} The scale of the polygon
 */
var Polygon = function(x, y, points, scale) {
    /**
     * The origin point of the polygon, all points are relative to this
     *
     * @property position
     * @type Vector
     */
    this.position = new Vector();

    /**
     * The unscaled points of the polygon, the X & Y values here should be
     * relative to the origin X & Y.
     *
     * @property _points
     * @type Array<Vector>
     * @default []
     * @private
     */
    this._points = null;

    /**
     * The scale of the polygon
     *
     * @property scale
     * @type Vector
     * @default new Vector(1, 1)
     */
    this.scale = scale || new Vector(1, 1);

    /**
     * The points of the polygon, the X & Y values here should be
     * relative to the origin X & Y values.
     *
     * @property points
     * @type Array<Vector>
     * @default []
     */
    this.points = [];

    /**
     * These vectors are calculated by `this.recalc()` and represent the edges
     * of the polygon defined by it's points.
     *
     * @property edges
     * @type Array<Vector>
     * @default []
     * @readOnly
     */
    this.edges = [];

    /**
     * These vectors are calculated by `this.recalc()` and represent the normals
     * of the polygon edges defined by it's points.
     *
     * @property normals
     * @type Array<Vector>
     * @default []
     * @readOnly
     */
    this.normals = [];

    //if this is a flat array of numbers, convert it to points
    if(typeof points[0] === 'number') {
        var p = [];
        for(var i = 0, il = points.length; i < il; i+=2) {
            p.push(
                new Vector(points[i], points[i + 1])
            );
        }

        points = p;
    }

    //assign the points
    this._points = points;

    //set position
    this.x = x || 0;
    this.y = y || 0;

    //recalculate scaled points, edges, and normals
    this.recalc();

    //internal shape type
    this._shapetype = C.SHAPE.POLYGON;
};

inherit(Polygon, Object, {
    /**
     * Creates a clone of this polygon
     *
     * @method clone
     * @return {Polygon} a copy of the polygon
     */
    clone: function() {
        var points = [];
        for (var i=0; i<this._points.length; i++) {
            points.push(this._points[i].clone());
        }

        return new Polygon(points, this.scale);
    },

    /**
     * Copies the values from another polygon to this one
     *
     * @method copy
     * @param polygon {Polygon} The polygon to copy vlaues from
     * @return {Polygon} Returns itself.
     * @chainable
     */
    copy: function(poly) {
        //copy the position
        this.position.copy(poly.position);

        //clone the points to this polygon
        this._points.length = this.points.length = 0;
        for(var i = 0; i < poly._points.length; ++i) {
            this._points.push(poly._points[i].clone());
        }

        this.scale.copy(poly.scale);

        //update our edges and normals
        this.recalc();

        return this;
    },

    /**
     * Checks if the x, and y coords passed to this function are contained within this polygon
     *
     * @method contains
     * @param x {Number} The X coord of the point to test
     * @param y {Number} The Y coord of the point to test
     * @return {Boolean} if the x/y coords are within this polygon
     */
    contains: function(x, y) {
        var inside = false;

        // use some raycasting to test hits
        // https://github.com/substack/point-in-polygon/blob/master/index.js
        for(var i = 0, j = this.points.length - 1; i < this.points.length; j = i++) {
            var xi = this.points[i].x, yi = this.points[i].y,
                xj = this.points[j].x, yj = this.points[j].y,
                intersect = ((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);

            if(intersect) inside = !inside;
        }

        return inside;
    },

    /**
     * Checks if this polygon's values are equal to anothers
     *
     * @method equals
     * @param polygon {Polygon} The polygon to check against
     * @return {Boolean} True if they are equal
     */
    equals: function(poly) {
        //check position and points array length
        if(!this.position.equals(poly.position) || this.points.length !== poly.points.length) {
            return false;
        }

        //check each point
        for(var i = 0; i < poly.points.length; ++i) {
            if(!this.points[i].equals(poly.points[i])) {
                return false;
            }
        }

        return true;
    },

    /**
     * Recalculates the scaled points, edges, and normals of this polygon
     * based on the relative points
     *
     * @method recalc
     * @return {Polygon} Returns itself.
     * @chainable
     */
    recalc: function() {
        var points = this._points,
            len = points.length,
            p1, p2, e, n, i = 0;

        //scale our points
        for(i = 0; i < len; i++) {
            if(!this.points[i])
                this.points[i] = new Vector();

            this.points[i].set(
                this._points[i].x * this.scale.x,
                this._points[i].y * this.scale.y
            );
        }

        // reset edges and normals
        this.edges.length = this.normals.length = 0;

        //calculate edges and normals
        for(i = 0; i < len; ++i) {
            p1 = points[i];
            p2 = i < len - 1 ? points[i + 1] : points[0];
            e = p2.clone().sub(p1);
            n = e.clone().perp().normalize();

            this.edges.push(e);
            this.normals.push(n);
        }

        return this;
    }
});

/**
 * The origin X coord of the polygon
 *
 * @property x
 * @type Number
 * @default 0
 */
Object.defineProperty(Polygon.prototype, 'x', {
    get: function() {
        return this.position.x;
    },
    set: function(v) {
        this.position.x = v;
    }
});

/**
 * The origin Y coord of the polygon
 *
 * @property x
 * @type Number
 * @default 0
 */
Object.defineProperty(Polygon.prototype, 'y', {
    get: function() {
        return this.position.y;
    },
    set: function(v) {
        this.position.y = v;
    }
});

module.exports = Polygon;

},{"../constants":11,"../math/Vector":47,"../utils/inherit":69}],37:[function(_dereq_,module,exports){
//var Rectangle = module.exports = require('pixi.js').Rectangle;

var inherit = _dereq_('../utils/inherit'),
    Polygon = _dereq_('./Polygon'),
    Vector = _dereq_('../math/Vector'),
    math = _dereq_('../math/math'),
    C = _dereq_('../constants');

/**
 * The Rectangle object is an area defined by its position, as indicated by its
 * top-left corner point (x, y) and by its width and its height.
 *
 * @class Rectangle
 * @constructor
 * @param x {Number} The X coord of the upper-left corner of the rectangle
 * @param y {Number} The Y coord of the upper-left corner of the rectangle
 * @param width {Number} The overall wisth of this rectangle
 * @param height {Number} The overall height of this rectangle
 */
var Rectangle = function(x, y, width, height) {
    /**
     * @property position
     * @type Vector
     * @default 0
     */
    this.position = new Vector();

    //set positon
    this.x = x || 0;
    this.y = y || 0;

    /**
     * @property _width
     * @type Number
     * @default 0
     * @private
     */
    this._width = width || 0;

    /**
     * @property _height
     * @type Number
     * @default 0
     * @private
     */
    this._height = height || 0;

    /**
     * @property halfWidth
     * @type Number
     * @default 0
     */
    this.halfWidth = this._width / 2;

    /**
     * @property halfHeight
     * @type Number
     * @default 0
     */
    this.halfHeight = this._height / 2;

    //internal shape type
    this._shapetype = C.SHAPE.RECTANGLE;
};

inherit(Rectangle, Object, {
    /**
     * Creates a clone of this Rectangle
     *
     * @method clone
     * @return {Rectangle} a copy of the rectangle
     */
    clone: function() {
        return new Rectangle(this.x, this.y, this._width, this._height);
    },

    /**
     * Copies the values from another rectangle to this one
     *
     * @method copy
     * @param rectangle {Rectangle} The rectangle to copy vlaues from
     * @return {Rectangle} Returns itself.
     * @chainable
     */
    copy: function(rect) {
        this.x = rect.x;
        this.y = rect.y;
        this.width = rect.width;
        this.height = rect.height;

        return this;
    },

    /**
     * Checks if the x, and y coords passed to this function are contained within this Rectangle
     *
     * @method contains
     * @param x {Number} The X coord of the point to test
     * @param y {Number} The Y coord of the point to test
     * @return {Boolean} if the x/y coords are within this Rectangle
     */
    contains: function(x, y) {
        if(this._width <= 0 || this._height <= 0)
            return false;

        var x1 = this.x;
        if(x >= x1 && x <= x1 + this._width) {
            var y1 = this.y;

            if(y >= y1 && y <= y1 + this._height) {
                return true;
            }
        }

        return false;
    },

    /**
     * Checks if this rectangle overlaps another
     *
     * @method overlaps
     * @param rect {Rectangle} The rectangle to check if this overlaps
     * @return {Boolean} if the rectangle overlaps
     */
    overlaps: function(rect) {
        return this.right > rect.x &&
                this.x < rect.right &&
                this.bottom > rect.y &&
                this.y < rect.bottom;
    },

    /**
     * Returns a polygon from this rectangle's points
     *
     * @method toPolygon
     * @return {Polygon} The new polygon
     */
    toPolygon: function(pos) {
        pos = pos || this.position;

        return new Polygon(this.x - pos.x, this.y - pos.y, [
            new Vector(pos.x, pos.y), //top-left
            new Vector(this.width, pos.y), //top-right
            new Vector(this.width, this.height), //bottom-right
            new Vector(pos.x, this.height) //bottom-left
        ]);
    },

    /**
     * Checks if this rectangle's values are equal to anothers
     *
     * @method equals
     * @param rectangle {Rectangle} The rectangle to check against
     * @return {Boolean} True if they are equal
     */
    equals: function(rect) {
        return this.position.equals(rect.position) &&
                this._width === rect._width &&
                this._height === rect._height;
    },

    /**
     * Combines two rectangles together to create a new rectangle
     *
     * @method union
     * @param rectangle {Rectangle} The rectangle to union with
     * @param [output] {Rectangle} The rectangle object to output to, a new one is created by default
     * @return {Rectangle} a new rectangle object that is the combonation of both
     */
    union: function(rect, out) {
        out = out || new Rectangle();

        out.x = math.min(this.x, rect.x);
        out.y = math.min(this.y, rect.y);
        out.width = math.max(this.right, rect.right);
        out.height = math.max(this.bottom, rect.bottom);

        return out;
    }
});

/**
 * The top-left X coord of the rectangle
 *
 * @property x
 * @type Number
 * @default 0
 */
Object.defineProperty(Rectangle.prototype, 'x', {
    get: function() {
        return this.position.x;
    },
    set: function(v) {
        this.position.x = v;
    }
});

/**
 * The top-left Y coord of the rectangle
 *
 * @property y
 * @type Number
 * @default 0
 */
Object.defineProperty(Rectangle.prototype, 'y', {
    get: function() {
        return this.position.y;
    },
    set: function(v) {
        this.position.y = v;
    }
});


/**
 * The width of the object
 *
 * @property width
 * @type Number
 * @default 0
 */
Object.defineProperty(Rectangle.prototype, 'width', {
    get: function() {
        return this._width;
    },
    set: function(w) {
        this._width = w || 0;
        this.halfWidth = this._width / 2;
    }
});

/**
 * The height of the object
 *
 * @property height
 * @type Number
 * @defualt 0
 */
Object.defineProperty(Rectangle.prototype, 'height', {
    get: function() {
        return this._height;
    },
    set: function(h) {
        this._height = h || 0;
        this.halfHeight = this._height / 2;
    }
});

/**
 * Returns the right most X coord
 *
 * @property right
 * @type Number
 * @readOnly
 */
Object.defineProperty(Rectangle.prototype, 'right', {
    get: function() {
        return this.x + this._width;
    }
});

/**
 * Returns the left most X coord
 *
 * @property left
 * @type Number
 * @readOnly
 */
Object.defineProperty(Rectangle.prototype, 'left', {
    get: function() {
        return this.x;
    }
});

/**
 * Returns the top most Y coord
 *
 * @property top
 * @type Number
 * @readOnly
 */
Object.defineProperty(Rectangle.prototype, 'top', {
    get: function() {
        return this.y;
    }
});

/**
 * Returns the bottom most Y coord
 *
 * @property bottom
 * @type Number
 * @readOnly
 */
Object.defineProperty(Rectangle.prototype, 'bottom', {
    get: function() {
        return this.y + this._height;
    }
});

/**
 * The perimeter of the rectangle
 *
 * @property perimeter
 * @type Number
 * @readOnly
 */
Object.defineProperty(Rectangle.prototype, 'perimeter', {
    get: function() {
        return 2 * (this._width + this._height);
    }
});

/**
 * The area of the rectangle
 *
 * @property area
 * @type Number
 * @readOnly
 */
Object.defineProperty(Rectangle.prototype, 'area', {
    get: function() {
        return this._width * this._height;
    }
});

module.exports = Rectangle;

},{"../constants":11,"../math/Vector":47,"../math/math":48,"../utils/inherit":69,"./Polygon":36}],38:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    Input = _dereq_('./Input'),
    GamepadButtons = _dereq_('./gamepad/GamepadButtons'),
    GamepadSticks = _dereq_('./gamepad/GamepadSticks');

/**
 * Controls input from gamepads
 *
 * @class Gamepad
 * @extends Input
 * @constructor
 */
var Gamepad = function() {
    Input.call(this);

    /**
     * Tracks if we are polling for status/connections
     *
     * @property ticking
     * @type Boolean
     * @readOnly
     */
    this.ticking = false;

    /**
     * The currently activated gamepads list
     *
     * @property pads
     * @type Array<Gamepad>
     * @readOnly
     */
    this.pads = [];

    /**
     * Timestamp tracking for state changes
     *
     * @property prevTimestamps
     * @type Array<Number>
     * @private
     */
    this.prevTimestamps = [];

    /**
     * Holds the button handler for gamepad button events
     *
     * @property buttons
     * @type GamepadButtons
     * @readOnly
     */
    this.buttons = new GamepadButtons();

    /**
     * Holds the stick handler for gamepad stick events
     *
     * @property sticks
     * @type GamepadSticks
     * @readOnly
     */
    this.sticks = new GamepadSticks();

    //Firefox uses connect/disconnect events so listen to those
    window.addEventListener('MozGamepadConnected', this.onGamepadConnect.bind(this), false);
    window.addEventListener('MozGamepadDisconnected', this.onGamepadDisconnect.bind(this), false);

    //Since chrome only supports polling, we have to start looping immediately
    if (!!navigator.webkitGamepads || !!navigator.webkitGetGamepads) {
        this.startPolling();
    }
};

inherit(Gamepad, Input, {
    /**
     * Called when a gamepad connects (FF Only)
     *
     * @method onGamepadDisconnect
     * @param event {GamepadConnectEvent}
     * @private
     */
    onGamepadConnect: function(event) {
        //add the gamepad to our list
        this.pads.push(event.gamepad);

        //start polling
        this.startPolling();
    },
    /**
     * Called when a gamepad disconnects (FF Only)
     *
     * @method onGamepadDisconnect
     * @param event {GamepadDisconnectEvent}
     * @private
     */
    onGamepadDisconnect: function(event) {
        //remove the gamepad from our list
        for(var i = 0, il = this.pads.length; i < il; ++i) {
            if(this.pads[i].index === event.gamepad.index) {
                this.pads.splice(i, 1);
                break;
            }
        }

        //if no pads left connected, stop polling
        if(this.pads.length === 0)
            this.stopPolling();
    },
    /**
     * Stats polling for new gamepads and status updates
     *
     * @method startPolling
     */
    startPolling: function() {
        if(this.ticking) return;

        this.ticking = true;
        this.update();
    },
    /**
     * Stops polling for new gamepads and status updates
     *
     * @method stopPolling
     */
    stopPolling: function() {
        this.ticking = false;
    },
    /**
     * Polls for newly connected gamepads (Chrome Only)
     *
     * @method pollGamepads
     */
    //called on Chrome, which doesn't do the connect/disconnect events
    pollGamepads: function() {
        //get a list of connected gamepads
        var rawPads = (navigator.webkitGetGamepads && navigator.webkitGetGamepads()) || navigator.webkitGamepads;

        if(rawPads) {
            //reset the pads list
            this.pads.length = 0;

            //don't use the raw array from the browser, since it can have "holes"
            //if you plug in 2, then remove the first the second one is still index 1 (not 0)
            for(var i = 0, il = rawPads.length; i < il; ++i) {
                if(rawPads[i]) {
                    this.pads.push(rawPads[i]);
                }
            }
        }
    },
    /**
     * Polls the gamepad object for status updates and emits events when they occur
     *
     * @method pollStatus
     */
    pollStatus: function() {
        for(var i = 0, il = this.pads.length; i < il; ++i) {
            var pad = this.pads[i];
            //don't do anything if the current timestamp is the same as the previous one
            //(meaning the state has not changed). This is a chrome-only feature right now,
            //so first we have to check if it is empty as well
            if(pad.timestamp && (pad.timestamp === this.prevTimestamps[i]))
                continue;

            this.prevTimestamps[i] = pad.timestamp;
            this.buttons.pollStatus(pad);
            this.sticks.pollStatus(pad);
        }
    },
    /**
     * Called each frame to update polling mechanisms
     *
     * @method update
     */
    update: function() {
        if(!this.ticking) return;

        //pollin' fo' pads
        this.pollGamepads();

        //poll for the status of our gamepads
        this.pollStatus();
    }
});

module.exports = Gamepad;

},{"../utils/inherit":69,"./Input":39,"./gamepad/GamepadButtons":43,"./gamepad/GamepadSticks":44}],39:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    EventEmitter = _dereq_('../utils/EventEmitter');

/**
 * The base Input object, holds common functions and properties between input types
 *
 * @class Input
 * @extends Object
 * @uses EventEmitter
 * @constructor
 * @param game {Game} The game instance this input belongs to
 */
var Input = function(game) {
    EventEmitter.call(this);

    /**
     * The game instance this input belongs to
     *
     * @property game
     * @type Game
     */
    this.game = game;
};

inherit(Input, Object, {
});

module.exports = Input;

},{"../utils/EventEmitter":64,"../utils/inherit":69}],40:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    Keyboard = _dereq_('./Keyboard'),
    Gamepad = _dereq_('./Gamepad'),
    Pointers = _dereq_('./Pointers');

/**
 * Manages all input handlers in a unified way
 *
 * @class InputManager
 * @extends Object
 * @constructor
 * @param game {Game} The game instance this input belongs to
 */
var InputManager = function(game) {
    /**
     * The game instance this manager belongs to
     *
     * @property game
     * @type Game
     */
    this.game = game;

    /**
     * The dom element to bind events to
     *
     * @property canvas
     * @type Game
     */
    this.canvas = game.canvas;

    /**
     * Holds the keyboard handler for keyboard events
     *
     * @property keyboard
     * @type Keyboard
     * @readOnly
     */
    this.keyboard = new Keyboard(game);

    /**
     * Holds the pointer handler for pointer events
     *
     * @property pointer
     * @type Pointer
     * @readOnly
     */
    this.pointers = new Pointers(game);

    /**
     * Holds the gamepad handler for gamepad events
     *
     * @property gamepad
     * @type Keyboard
     * @readOnly
     */
    this.gamepad = new Gamepad();
};

inherit(InputManager, Object, {
    /**
     * Called internally every frame. Updates all the pointers and gamepad
     *
     * @method update
     * @param dt {Number} The delta time (in seconds) since the last update
     * @private
     */
    update: function(dt) {
        this.pointers.update(dt);
        this.gamepad.update(dt);
    }
});

module.exports = InputManager;

},{"../utils/inherit":69,"./Gamepad":38,"./Keyboard":41,"./Pointers":42}],41:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    Input = _dereq_('./Input');

/**
 * Controls keyboard input
 *
 * @class Keyboard
 * @extends Input
 * @constructor
 * @param game {Game} The game instance this input belongs to
 */
var Keyboard = function(game) {
    Input.call(this, game);

    /**
     * The current sequence of keys that have been pressed
     *
     * @property sequence
     * @type Array<Number>
     * @readOnly
     */
    this.sequence = [];

    /**
     * The amount of time it takes for the sequence to clear out, in ms
     *
     * @property sequenceTimeout
     * @type Number
     * @default 500
     */
    this.sequenceTimeout = 500;

    /**
     * The timeout ID for the wait to clear the input sequence
     *
     * @property _clearSq
     * @type Number
     * @private 
     */
    this._clearSq = null;

    game.canvas.addEventListener('keydown', this.onKeyDown.bind(this), false);
    game.canvas.addEventListener('keyup', this.onKeyUp.bind(this), false);
};

inherit(Keyboard, Input, {
    /**
     * Called when a key is pressed down
     *
     * @method onKeyDown
     * @param event {DOMEvent}
     * @param override {Number} The key code to use instead of checking event data
     * @private
     */
    onKeyDown: function(e, override) {
        //if(e.target === this.view.parentElement)
        return this.modifyKey(e, override || e.keyCode || e.which, true);
    },
    /**
     * Called when a key is released
     *
     * @method onKeyUp
     * @param event {DOMEvent}
     * @param override {Number} The key code to use instead of checking event data
     * @private
     */
    onKeyUp: function(e, override) {
        //if(e.target === this.view.parentElement)
        return this.modifyKey(e, override || e.keyCode || e.which, false);
    },
    /**
     * Called when a key state has changed, updates current sequence and emits events
     *
     * @method modifyKey
     * @param event {DOMEvent}
     * @param key {Number} The key code that has changed
     * @param down {Boolean} Whether the key has been pressed or not
     * @private
     */
    modifyKey: function(e, key, down) {
        //emit single key event
        this.emit(key, this._getEventData(e, down));

        //when pressed is when we process a key for a sequence
        if(down) {
            //update the key sequence
            this.sequence.push(key);

            //process current sequence
            var s = this.sequence.toString();
            if(s !== key.toString()) {
                this.emit(s, this._getEventData(e, down));
            }

            //set timeout to clear sequence
            clearTimeout(this._clearSq);
            this._clearSq = setTimeout(this._clearSequence.bind(this), this.sequenceTimeout);
        }
    },
    /**
     * Generates an event data object for a keyboard event
     *
     * @method _getEventData
     * @param event {DOMEvent} The original DOMEvent that was passed into the raw event handler
     * @param down {Boolean} Is this a keydown event
     * @return {Object} The event object
     * @private
     */
    _getEventData: function(e, down) {
        return {
            input: this,
            originalEvent: e,
            down: down
        };
    },
    /**
     * Clears the current sequence so that a new one can start
     *
     * @method _clearSequence
     * @private
     */
    _clearSequence: function() {
        this.sequence.length = 0;
    }
});

/**
 * Bindable keycodes
 *
 * @property KEY
 * @type Object
 * @static
 */
Keyboard.KEY = {
    BACKSPACE: 8,
    TAB: 9,
    ENTER: 13,
    SHIFT: 16,
    CTRL: 17,
    ALT: 18,
    PAUSE: 19,
    ESC: 27,
    SPACE: 32,
    PAGE_UP: 33,
    PAGE_DOWN: 34,
    END: 35,
    HOME: 36,
    LEFT: 37,
    UP: 38,
    RIGHT: 39,
    DOWN: 40,
    INSERT: 45,
    DELETE: 46,
    NUM0: 48,
    NUM1: 49,
    NUM2: 50,
    NUM3: 51,
    NUM4: 52,
    NUM5: 53,
    NUM6: 54,
    NUM7: 55,
    NUM8: 56,
    NUM9: 57,
    PLUS: 61,
    A: 65,
    B: 66,
    C: 67,
    D: 68,
    E: 69,
    F: 70,
    G: 71,
    H: 72,
    I: 73,
    J: 74,
    K: 75,
    L: 76,
    M: 77,
    N: 78,
    O: 79,
    P: 80,
    Q: 81,
    R: 82,
    S: 83,
    T: 84,
    U: 85,
    V: 86,
    W: 87,
    X: 88,
    Y: 89,
    Z: 90,
    NUMPAD0: 96,
    NUMPAD1: 97,
    NUMPAD2: 98,
    NUMPAD3: 99,
    NUMPAD4: 100,
    NUMPAD5: 101,
    NUMPAD6: 102,
    NUMPAD7: 103,
    NUMPAD8: 104,
    NUMPAD9: 105,
    NUMPAD_STAR: 106,
    NUMPAD_PLUS: 107,
    NUMPAD_MINUS: 109,
    NUMPAD_DOT: 110,
    NUMPAD_SLASH: 111,
    F1: 112,
    F2: 113,
    F3: 114,
    F4: 115,
    MINUS: 173,
    TILDE: 192
};

module.exports = Keyboard;

},{"../utils/inherit":69,"./Input":39}],42:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    Input = _dereq_('./Input'),
    Pointer = _dereq_('./pointer/Pointer');

/**
 * Controls pointer input (mouse, touch, pen, etc) or all pointers tracked by the game
 *
 * @class Pointers
 * @extends Input
 * @constructor
 * @param game {Game} The game instance this input belongs to
 */
//TODO: Sprite interactivity and Interaction History
var Pointers = function(game) {
    Input.call(this, game);

    /**
     * The pointer instances currently being used, keyed by an ID
     *
     * @property pointers
     * @type Object<Pointer>
     */
    this.pointers = {};

    /**
     * The max number of pointers to track
     *
     * @property maxPointers
     * @type Number
     * @default 10
     */
    this.maxPointers = 10;

    /**
     * The time that must pass between a down (touchstart/mousedown) and up (touchend/mouseup)
     * event for it to be considered a "click" event, in milliseconds
     *
     * @property clickDelay
     * @type Number
     * @default 200
     */
    this.clickDelay = 200;

    /**
     * The max time that can pass between two click events for it to be considered a
     * "doubleclick" event, in milliseconds
     *
     * @property doubleClickDelay
     * @type Number
     * @default 300
     */
    this.doubleClickDelay = 300;

    /**
     * The time that must pass after a down event for it to be considered a "hold" event, in milliseconds
     *
     * @property holdDelay
     * @type Number
     * @default 2000
     */
    this.holdDelay = 2000;

    //create the mouse pointer object
    this.mouse = this.pointers[1] = new Pointer(1, this);

    //number of pointers being tracked
    this.activePointers = 0;

    //bind events
    game.canvas.addEventListener('pointerdown',     this.onPointer.bind(this, 'down'),    false);
    game.canvas.addEventListener('pointerup',       this.onPointer.bind(this, 'up'),      false);
    game.canvas.addEventListener('pointermove',     this.onPointer.bind(this, 'move'),    false);
    game.canvas.addEventListener('pointerover',     this.onPointer.bind(this, 'over'),    false);
    game.canvas.addEventListener('pointerout',      this.onPointer.bind(this, 'out'),     false);
    game.canvas.addEventListener('pointercancel',   this.onPointer.bind(this, 'cancel'),  false);
    game.canvas.addEventListener('pointerenter',    this.onPointer.bind(this, 'enter'),   false);
    game.canvas.addEventListener('pointerleave',    this.onPointer.bind(this, 'leave'),   false);

    /**
     * Fired when a pointer is pressed on the canvas
     *
     * @event down
     * @param pointer {Pointer} The pointer instance that had a 'pointerdown' event
     */

    /**
     * Fired when a pointer is released off the canvas
     *
     * @event up
     * @param pointer {Pointer} The pointer instance that had a 'pointerup' event
     */

    /**
     * Fired when a pointer is moved while on the canvas
     *
     * @event move
     * @param pointer {Pointer} The pointer instance that had a 'pointermove' event
     */

    /**
     * Fired when a pointer moves over the canvas
     *
     * @event over
     * @param pointer {Pointer} The pointer instance that had a 'pointerover' event
     */

    /**
     * Fired when a pointer moves out of the canvas
     *
     * @event out
     * @param pointer {Pointer} The pointer instance that had a 'pointerout' event
     */

    /**
     * Fired when a pointer event is canceled
     *
     * @event cancel
     * @param pointer {Pointer} The pointer instance that had a 'pointercancel' event
     */

    /**
     * Fired when a pointer enters the canvas
     *
     * @event enter
     * @param pointer {Pointer} The pointer instance that had a 'pointerenter' event
     */

    /**
     * Fired when a pointer leaves the canvas
     *
     * @event leave
     * @param pointer {Pointer} The pointer instance that had a 'pointerleave' event
     */

    this.interactiveSprites = [];
};

inherit(Pointers, Input, {
    watchSprite: function(spr) {
        this.interactiveSprites.push(spr);
        return this;
    },
    /**
     * Callback that is called when a pointer event occurs.
     *
     * @method onPointer
     * @param name {String} The name of the pointer event with out the 'pointer' prefix
     * @param evt {DOMEvent} The DOM Event
     * @private
     */
    onPointer: function(name, evt) {
        var id = evt.pointerId,
            pointer = this.pointers[id];

        //create a new pointer object if we need it, and if there is room
        //if there isn't room for a new pointer object then we just return
        //without echoing the event or anything.
        if(!pointer) {
            if(this._numPointers < this.maxPointers) {
                this.pointers[id] = new Pointer(id, this);
            } else {
                return;
            }
        }

        if(pointer[name])
            pointer[name](evt);

        this.emit(name, pointer);
    },
    /**
     * Called internally every frame. Updates all the pointers
     *
     * @method update
     * @param dt {Number} The delta time (in seconds) since the last update
     * @return {Pointers} Returns iteself for chainability
     * @private
     */
    update: function(dt) {
        var pointer;

        //go through each pointer and update
        for(var id in this.pointers) {
            pointer = this.pointers[id];

            if(pointer) {
                pointer.update(dt);
            }
        }

        return this;
    }
});

module.exports = Pointers;

//////////////////////////////////////////////////////////////////////////////////////////////////
// hand-1.1.2.js
//////////////////////////////////////////////////////////////////////////////////////////////////

/* jshint ignore:start */

(function () {
    // Polyfilling indexOf for old browsers
    if (!Array.prototype.indexOf) {
        Array.prototype.indexOf = function (searchElement) {
            var t = Object(this);
            var len = t.length >>> 0;
            if (len === 0) {
                return -1;
            }
            var n = 0;
            if (arguments.length > 0) {
                n = Number(arguments[1]);
                if (n != n) { // shortcut for verifying if it's NaN
                    n = 0;
                } else if (n != 0 && n != Infinity && n != -Infinity) {
                    n = (n > 0 || -1) * Math.floor(Math.abs(n));
                }
            }
            if (n >= len) {
                return -1;
            }
            var k = n >= 0 ? n : Math.max(len - Math.abs(n), 0);
            for (; k < len; k++) {
                if (k in t && t[k] === searchElement) {
                    return k;
                }
            }
            return -1;
        };
    }

    // Installing Hand.js
    var supportedEventsNames = ["PointerDown", "PointerUp", "PointerMove", "PointerOver", "PointerOut", "PointerCancel", "PointerEnter", "PointerLeave",
                                "pointerdown", "pointerup", "pointermove", "pointerover", "pointerout", "pointercancel", "pointerenter", "pointerleave"
    ];

    var POINTER_TYPE_TOUCH = "touch";
    var POINTER_TYPE_PEN = "pen";
    var POINTER_TYPE_MOUSE = "mouse";

    var previousTargets = {};

    var checkPreventDefault = function(element) {
        while (element && element.handjs_forcePreventDefault !== true) {
            element = element.parentElement;
        }
        return element != null;
    };

    // Touch events
    var generateTouchClonedEvent = function (sourceEvent, newName) {
        // Considering touch events are almost like super mouse events
        var evObj;

        if (document.createEvent) {
            evObj = document.createEvent('MouseEvents');
            evObj.initMouseEvent(newName, true, true, window, 1, sourceEvent.screenX, sourceEvent.screenY,
                sourceEvent.clientX, sourceEvent.clientY, sourceEvent.ctrlKey, sourceEvent.altKey,
                sourceEvent.shiftKey, sourceEvent.metaKey, sourceEvent.button, null);
        }
        else {
            evObj = document.createEventObject();
            evObj.screenX = sourceEvent.screenX;
            evObj.screenY = sourceEvent.screenY;
            evObj.clientX = sourceEvent.clientX;
            evObj.clientY = sourceEvent.clientY;
            evObj.ctrlKey = sourceEvent.ctrlKey;
            evObj.altKey = sourceEvent.altKey;
            evObj.shiftKey = sourceEvent.shiftKey;
            evObj.metaKey = sourceEvent.metaKey;
            evObj.button = sourceEvent.button;
        }
        // offsets
        if (evObj.offsetX === undefined) {
            if (sourceEvent.offsetX !== undefined) {

                // For Opera which creates readonly properties
                if (Object && Object.defineProperty !== undefined) {
                    Object.defineProperty(evObj, "offsetX", {
                        writable: true
                    });
                    Object.defineProperty(evObj, "offsetY", {
                        writable: true
                    });
                }

                evObj.offsetX = sourceEvent.offsetX;
                evObj.offsetY = sourceEvent.offsetY;
            }
            else if (sourceEvent.layerX !== undefined) {
                evObj.offsetX = sourceEvent.layerX - sourceEvent.currentTarget.offsetLeft;
                evObj.offsetY = sourceEvent.layerY - sourceEvent.currentTarget.offsetTop;
            }
        }

        // adding missing properties

        if (sourceEvent.isPrimary !== undefined)
            evObj.isPrimary = sourceEvent.isPrimary;
        else
            evObj.isPrimary = true;

        if (sourceEvent.pressure)
            evObj.pressure = sourceEvent.pressure;
        else {
            var button = 0;

            if (sourceEvent.which !== undefined)
                button = sourceEvent.which;
            else if (sourceEvent.button !== undefined) {
                button = sourceEvent.button;
            }
            evObj.pressure = (button == 0) ? 0 : 0.5;
        }


        if (sourceEvent.rotation)
            evObj.rotation = sourceEvent.rotation;
        else
            evObj.rotation = 0;

        // Timestamp
        if (sourceEvent.hwTimestamp)
            evObj.hwTimestamp = sourceEvent.hwTimestamp;
        else
            evObj.hwTimestamp = 0;

        // Tilts
        if (sourceEvent.tiltX)
            evObj.tiltX = sourceEvent.tiltX;
        else
            evObj.tiltX = 0;

        if (sourceEvent.tiltY)
            evObj.tiltY = sourceEvent.tiltY;
        else
            evObj.tiltY = 0;

        // Width and Height
        if (sourceEvent.height)
            evObj.height = sourceEvent.height;
        else
            evObj.height = 0;

        if (sourceEvent.width)
            evObj.width = sourceEvent.width;
        else
            evObj.width = 0;

        // preventDefault
        evObj.preventDefault = function () {
            if (sourceEvent.preventDefault !== undefined)
                sourceEvent.preventDefault();
        };

        // stopPropagation
        if (evObj.stopPropagation !== undefined) {
            var current = evObj.stopPropagation;
            evObj.stopPropagation = function () {
                if (sourceEvent.stopPropagation !== undefined)
                    sourceEvent.stopPropagation();
                current.call(this);
            };
        }

        // Constants
        evObj.POINTER_TYPE_TOUCH = POINTER_TYPE_TOUCH;
        evObj.POINTER_TYPE_PEN = POINTER_TYPE_PEN;
        evObj.POINTER_TYPE_MOUSE = POINTER_TYPE_MOUSE;

        // Pointer values
        evObj.pointerId = sourceEvent.pointerId;
        evObj.pointerType = sourceEvent.pointerType;

        switch (evObj.pointerType) {// Old spec version check
            case 2:
                evObj.pointerType = evObj.POINTER_TYPE_TOUCH;
                break;
            case 3:
                evObj.pointerType = evObj.POINTER_TYPE_PEN;
                break;
            case 4:
                evObj.pointerType = evObj.POINTER_TYPE_MOUSE;
                break;
        }

        // If force preventDefault
        if (sourceEvent.currentTarget && checkPreventDefault(sourceEvent.currentTarget) === true) {
            evObj.preventDefault();
        }

        // Fire event
        if (sourceEvent.target) {
            sourceEvent.target.dispatchEvent(evObj);
        } else {
            sourceEvent.srcElement.fireEvent("on" + getMouseEquivalentEventName(newName), evObj); // We must fallback to mouse event for very old browsers
        }
    };

    var generateMouseProxy = function (evt, eventName) {
        evt.pointerId = 1;
        evt.pointerType = POINTER_TYPE_MOUSE;
        generateTouchClonedEvent(evt, eventName);
    };

    var generateTouchEventProxy = function (name, touchPoint, target, eventObject) {
        var touchPointId = touchPoint.identifier + 2; // Just to not override mouse id

        touchPoint.pointerId = touchPointId;
        touchPoint.pointerType = POINTER_TYPE_TOUCH;
        touchPoint.currentTarget = target;
        touchPoint.target = target;

        if (eventObject.preventDefault !== undefined) {
            touchPoint.preventDefault = function () {
                eventObject.preventDefault();
            };
        }

        generateTouchClonedEvent(touchPoint, name);
    };

    var checkRegisteredEvents = function(element, eventName) {
        while (element && !(element.__handjsGlobalRegisteredEvents && element.__handjsGlobalRegisteredEvents[eventName])) {
            element = element.parentElement;
        }
        return element != null;
    };

    var generateTouchEventProxyIfRegistered = function (eventName, touchPoint, target, eventObject) { // Check if user registered this event
        if (checkRegisteredEvents(target, eventName)) {
            generateTouchEventProxy(eventName, touchPoint, target, eventObject);
        }
    };

    var handleOtherEvent = function (eventObject, name, useLocalTarget, checkRegistration) {
        if (eventObject.preventManipulation)
            eventObject.preventManipulation();

        for (var i = 0; i < eventObject.changedTouches.length; ++i) {
            var touchPoint = eventObject.changedTouches[i];

            if (useLocalTarget) {
                previousTargets[touchPoint.identifier] = touchPoint.target;
            }

            if (checkRegistration) {
                generateTouchEventProxyIfRegistered(name, touchPoint, previousTargets[touchPoint.identifier], eventObject);
            } else {
                generateTouchEventProxy(name, touchPoint, previousTargets[touchPoint.identifier], eventObject);
            }
        }
    };

    var getMouseEquivalentEventName = function (eventName) {
        return eventName.toLowerCase().replace("pointer", "mouse");
    };

    var getPrefixEventName = function (item, prefix, eventName) {
        var newEventName;

        if (eventName == eventName.toLowerCase()) {
            var indexOfUpperCase = supportedEventsNames.indexOf(eventName) - (supportedEventsNames.length / 2);
            newEventName = prefix + supportedEventsNames[indexOfUpperCase];
        }
        else {
            newEventName = prefix + eventName;
        }

        // Fallback to PointerOver if PointerEnter is not currently supported
        if (newEventName === prefix + "PointerEnter" && item["on" + prefix.toLowerCase() + "pointerenter"] === undefined) {
            newEventName = prefix + "PointerOver";
        }

        // Fallback to PointerOut if PointerLeave is not currently supported
        if (newEventName === prefix + "PointerLeave" && item["on" + prefix.toLowerCase() + "pointerleave"] === undefined) {
            newEventName = prefix + "PointerOut";
        }

        return newEventName;
    };

    var registerOrUnregisterEvent = function (item, name, func, enable) {
        if (item.__handjsRegisteredEvents === undefined) {
            item.__handjsRegisteredEvents = [];
        }

        if (enable) {
            if (item.__handjsRegisteredEvents[name] !== undefined) {
                item.__handjsRegisteredEvents[name]++;
                return;
            }

            item.__handjsRegisteredEvents[name] = 1;
            item.addEventListener(name, func, false);
        } else {

            if (item.__handjsRegisteredEvents.indexOf(name) !== -1) {
                item.__handjsRegisteredEvents[name]--;

                if (item.__handjsRegisteredEvents[name] != 0) {
                    return;
                }
            }
            item.removeEventListener(name, func);
            item.__handjsRegisteredEvents[name] = 0;
        }
    };

    var setTouchAware = function (item, eventName, enable) {
        // If item is already touch aware, do nothing
        if (item.onpointerdown !== undefined) {
            return;
        }

        // IE 10
        if (item.onmspointerdown !== undefined) {
            var msEventName = getPrefixEventName(item, "MS", eventName);

            registerOrUnregisterEvent(item, msEventName, function (evt) { generateTouchClonedEvent(evt, eventName); }, enable);

            // We can return because MSPointerXXX integrate mouse support
            return;
        }

        // Chrome, Firefox
        if (item.ontouchstart !== undefined) {
            switch (eventName) {
                case "pointermove":
                    registerOrUnregisterEvent(item, "touchmove", function (evt) { handleOtherEvent(evt, eventName); }, enable);
                    break;
                case "pointercancel":
                    registerOrUnregisterEvent(item, "touchcancel", function (evt) { handleOtherEvent(evt, eventName); }, enable);
                    break;
                case "pointerdown":
                case "pointerup":
                case "pointerout":
                case "pointerover":
                case "pointerleave":
                case "pointerenter":
                    // These events will be handled by the window.ontouchmove function
                    if (!item.__handjsGlobalRegisteredEvents) {
                        item.__handjsGlobalRegisteredEvents = [];
                    }

                    if (enable) {
                        if (item.__handjsGlobalRegisteredEvents[eventName] !== undefined) {
                            item.__handjsGlobalRegisteredEvents[eventName]++;
                            return;
                        }
                        item.__handjsGlobalRegisteredEvents[eventName] = 1;
                    } else {
                        if (item.__handjsGlobalRegisteredEvents[eventName] !== undefined) {
                            item.__handjsGlobalRegisteredEvents[eventName]--;
                            if (item.__handjsGlobalRegisteredEvents[eventName] < 0) {
                                item.__handjsGlobalRegisteredEvents[eventName] = 0;
                            }
                        }
                    }
                    break;
            }
        }

        // Fallback to mouse
        switch (eventName) {
            case "pointerdown":
                registerOrUnregisterEvent(item, "mousedown", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                break;
            case "pointermove":
                registerOrUnregisterEvent(item, "mousemove", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                break;
            case "pointerup":
                registerOrUnregisterEvent(item, "mouseup", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                break;
            case "pointerover":
                registerOrUnregisterEvent(item, "mouseover", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                break;
            case "pointerout":
                registerOrUnregisterEvent(item, "mouseout", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                break;
            case "pointerenter":
                if (item.onmouseenter === undefined) { // Fallback to mouseover
                    registerOrUnregisterEvent(item, "mouseover", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                } else {
                    registerOrUnregisterEvent(item, "mouseenter", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                }
                break;
            case "pointerleave":
                if (item.onmouseleave === undefined) { // Fallback to mouseout
                    registerOrUnregisterEvent(item, "mouseout", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                } else {
                    registerOrUnregisterEvent(item, "mouseleave", function (evt) { generateMouseProxy(evt, eventName); }, enable);
                }
                break;
        }
    };

    // Intercept addEventListener calls by changing the prototype
    var interceptAddEventListener = function (root) {
        var current = root.prototype ? root.prototype.addEventListener : root.addEventListener;

        var customAddEventListener = function (name, func, capture) {
            // Branch when a PointerXXX is used
            if (supportedEventsNames.indexOf(name) != -1) {
                setTouchAware(this, name, true);
            }

            if (current === undefined) {
                this.attachEvent("on" + getMouseEquivalentEventName(name), func);
            } else {
                current.call(this, name, func, capture);
            }
        };

        if (root.prototype) {
            root.prototype.addEventListener = customAddEventListener;
        } else {
            root.addEventListener = customAddEventListener;
        }
    };

    // Intercept removeEventListener calls by changing the prototype
    var interceptRemoveEventListener = function (root) {
        var current = root.prototype ? root.prototype.removeEventListener : root.removeEventListener;

        var customRemoveEventListener = function (name, func, capture) {
            // Release when a PointerXXX is used
            if (supportedEventsNames.indexOf(name) != -1) {
                setTouchAware(this, name, false);
            }

            if (current === undefined) {
                this.detachEvent(getMouseEquivalentEventName(name), func);
            } else {
                current.call(this, name, func, capture);
            }
        };
        if (root.prototype) {
            root.prototype.removeEventListener = customRemoveEventListener;
        } else {
            root.removeEventListener = customRemoveEventListener;
        }
    };

    // Hooks
    interceptAddEventListener(HTMLElement);
    interceptAddEventListener(document);
    interceptAddEventListener(HTMLBodyElement);
    interceptAddEventListener(HTMLDivElement);
    interceptAddEventListener(HTMLImageElement);
    interceptAddEventListener(HTMLUListElement);
    interceptAddEventListener(HTMLAnchorElement);
    interceptAddEventListener(HTMLLIElement);
    interceptAddEventListener(HTMLTableElement);
    if (window.HTMLSpanElement) {
        interceptAddEventListener(HTMLSpanElement);
    }
    if (window.HTMLCanvasElement) {
        interceptAddEventListener(HTMLCanvasElement);
    }
    if (window.SVGElement) {
        interceptAddEventListener(SVGElement);
    }

    interceptRemoveEventListener(HTMLElement);
    interceptRemoveEventListener(document);
    interceptRemoveEventListener(HTMLBodyElement);
    interceptRemoveEventListener(HTMLDivElement);
    interceptRemoveEventListener(HTMLImageElement);
    interceptRemoveEventListener(HTMLUListElement);
    interceptRemoveEventListener(HTMLAnchorElement);
    interceptRemoveEventListener(HTMLLIElement);
    interceptRemoveEventListener(HTMLTableElement);
    if (window.HTMLSpanElement) {
        interceptRemoveEventListener(HTMLSpanElement);
    }
    if (window.HTMLCanvasElement) {
        interceptRemoveEventListener(HTMLCanvasElement);
    }
    if (window.SVGElement) {
        interceptRemoveEventListener(SVGElement);
    }

    // Handling move on window to detect pointerleave/out/over
    if (window.ontouchstart !== undefined) {
        window.addEventListener('touchstart', function (eventObject) {
            for (var i = 0; i < eventObject.changedTouches.length; ++i) {
                var touchPoint = eventObject.changedTouches[i];
                previousTargets[touchPoint.identifier] = touchPoint.target;

                generateTouchEventProxyIfRegistered("pointerenter", touchPoint, touchPoint.target, eventObject);
                generateTouchEventProxyIfRegistered("pointerover", touchPoint, touchPoint.target, eventObject);
                generateTouchEventProxyIfRegistered("pointerdown", touchPoint, touchPoint.target, eventObject);
            }
        });

        window.addEventListener('touchend', function (eventObject) {
            for (var i = 0; i < eventObject.changedTouches.length; ++i) {
                var touchPoint = eventObject.changedTouches[i];
                var currentTarget = previousTargets[touchPoint.identifier];

                generateTouchEventProxyIfRegistered("pointerup", touchPoint, currentTarget, eventObject);
                generateTouchEventProxyIfRegistered("pointerout", touchPoint, currentTarget, eventObject);
                generateTouchEventProxyIfRegistered("pointerleave", touchPoint, currentTarget, eventObject);
            }
        });

        window.addEventListener('touchmove', function (eventObject) {
            for (var i = 0; i < eventObject.changedTouches.length; ++i) {
                var touchPoint = eventObject.changedTouches[i];
                var newTarget = document.elementFromPoint(touchPoint.clientX, touchPoint.clientY);
                var currentTarget = previousTargets[touchPoint.identifier];

                if (currentTarget === newTarget) {
                    continue; // We can skip this as the pointer is effectively over the current target
                }

                if (currentTarget) {
                    // Raise out
                    generateTouchEventProxyIfRegistered("pointerout", touchPoint, currentTarget, eventObject);

                    // Raise leave
                    if (!currentTarget.contains(newTarget)) { // Leave must be called if the new target is not a child of the current
                        generateTouchEventProxyIfRegistered("pointerleave", touchPoint, currentTarget, eventObject);
                    }
                }

                if (newTarget) {
                    // Raise over
                    generateTouchEventProxyIfRegistered("pointerover", touchPoint, newTarget, eventObject);

                    // Raise enter
                    if (!newTarget.contains(currentTarget)) { // Leave must be called if the new target is not the parent of the current
                        generateTouchEventProxyIfRegistered("pointerenter", touchPoint, newTarget, eventObject);
                    }
                }
                previousTargets[touchPoint.identifier] = newTarget;
            }
        });
    }

    // Extension to navigator
    if (navigator.pointerEnabled === undefined) {

        // Indicates if the browser will fire pointer events for pointing input
        navigator.pointerEnabled = true;

        // IE
        if (navigator.msPointerEnabled) {
            navigator.maxTouchPoints = navigator.msMaxTouchPoints;
        }
    }

    // Handling touch-action css rule
    if (document.styleSheets && document.addEventListener) {
        document.addEventListener("DOMContentLoaded", function () {

            var trim = function (string) {
                return string.replace(/^\s+|\s+$/, '');
            };

            var processStylesheet = function (unfilteredSheet) {
                var globalRegex = new RegExp(".+?{.*?}", "m");
                var selectorRegex = new RegExp(".+?{", "m");

                while (unfilteredSheet != "") {
                    var filter = globalRegex.exec(unfilteredSheet);
                    if (!filter) {
                        break;
                    }
                    var block = filter[0];
                    unfilteredSheet = trim(unfilteredSheet.replace(block, ""));
                    var selectorText = trim(selectorRegex.exec(block)[0].replace("{", ""));

                    // Checking if the user wanted to deactivate the default behavior
                    if (block.replace(/\s/g, "").indexOf("touch-action:none") != -1) {
                        var elements = document.querySelectorAll(selectorText);

                        for (var elementIndex = 0; elementIndex < elements.length; elementIndex++) {
                            var element = elements[elementIndex];

                            if (element.style.msTouchAction !== undefined) {
                                element.style.msTouchAction = "none";
                            }
                            else {
                                element.handjs_forcePreventDefault = true;
                            }
                        }
                    }
                }
            }; // Looking for touch-action in referenced stylesheets
            try {
                for (var index = 0; index < document.styleSheets.length; index++) {
                    var sheet = document.styleSheets[index];

                    if (sheet.href == undefined) { // it is an inline style
                        continue;
                    }

                    // Loading the original stylesheet
                    var xhr = new XMLHttpRequest();
                    xhr.open("get", sheet.href, false);
                    xhr.send();

                    var unfilteredSheet = xhr.responseText.replace(/(\n|\r)/g, "");

                    processStylesheet(unfilteredSheet);
                }
            } catch (e) {
                // Silently fail...
            }

            // Looking for touch-action in inline styles
            var styles = document.getElementsByTagName("style");
            for (var index = 0; index < styles.length; index++) {
                var inlineSheet = styles[index];

                var inlineUnfilteredSheet = trim(inlineSheet.innerHTML.replace(/(\n|\r)/g, ""));

                processStylesheet(inlineUnfilteredSheet);
            }
        }, false);
    }

})();

/* jshint ignore:end */

},{"../utils/inherit":69,"./Input":39,"./pointer/Pointer":45}],43:[function(_dereq_,module,exports){
var inherit = _dereq_('../../utils/inherit'),
    Input = _dereq_('../Input');

/**
 * Controls gamepad button input
 *
 * @class GamepadButtons
 * @extends Input
 * @constructor
 */
var GamepadButtons = function() {
    Input.call(this);

    /**
     * The threshold at which we consider a button "pressed"
     *
     * @property threshold
     * @type Number
     * @default 0.4
     */
    this.threshold = 0.4;

    /**
     * Track the status of each button on the gamepad
     *
     * @property buttons
     * @type Object
     */
    this.buttons = {};

    //setup default objects for each axis
    for(var bt in GamepadButtons.BUTTON) {
        this.buttons[GamepadButtons.BUTTON[bt]] = {
            code: GamepadButtons.BUTTON[bt],
            name: bt,
            down: false,
            value: 0
        };
    }
};

inherit(GamepadButtons, Input, {
    /**
     * Polls the gamepad object for status updates and emits events when they occur
     *
     * @method pollStatus
     * @param pad {Gamepad} The gamepad object to check
     */
    pollStatus: function(pad) {
        for(var b = 0, bl = pad.buttons.length; b < bl; ++b) {
            var down = (pad.buttons[b] > this.threshold),
                status = this.buttons[b];

            status.value = pad.buttons[b];

            //down state changed
            if(status.down !== down) {
                status.down = down;

                this.emit(b, status);
            }
        }
    }
});

/**
 * Bindable Gamepad Buttons
 *
 * @property BUTTON
 * @type Object
 * @static
 */
GamepadButtons.BUTTON = {
    FACE_1: 0, // Face (main) buttons
    FACE_2: 1,
    FACE_3: 2,
    FACE_4: 3,
    LEFT_SHOULDER: 4, // Top shoulder buttons
    RIGHT_SHOULDER: 5,
    LEFT_TRIGGER: 6, // Bottom shoulder buttons
    RIGHT_TRIGGER: 7,
    SELECT: 8,
    START: 9,
    LEFT_ANALOGUE_STICK: 10, // Analogue sticks (if depressible)
    RIGHT_ANALOGUE_STICK: 11,
    PAD_UP: 12, // Directional (discrete) pad
    PAD_DOWN: 13,
    PAD_LEFT: 14,
    PAD_RIGHT: 15,
    SYSTEM_MENU: 16   // on console controllers this would be the button to open the system menu
};

GamepadButtons.getGpButtonName = function(i) {
    for(var k in GamepadButtons.BUTTON) {
        if(GamepadButtons.BUTTON[k] === i) {
            return k;
        }
    }

    return '';
};

module.exports = GamepadButtons;

},{"../../utils/inherit":69,"../Input":39}],44:[function(_dereq_,module,exports){
var inherit = _dereq_('../../utils/inherit'),
    Input = _dereq_('../Input');

/**
 * Controls gamepad stick input
 *
 * @class GamepadSticks
 * @extends Input
 * @constructor
 */
var GamepadSticks = function() {
    Input.call(this);

    /**
     * The threshold at which we consider a stick moved from center
     *
     * @property threshold
     * @type Number
     * @default 0.5
     */
    this.threshold = 0.5;

    /**
     * Track the status of each of the axes on the gamepad
     *
     * @property axes
     * @type Object
     */
    this.axes = {};

    //setup default objects for each axis
    for(var ax in GamepadSticks.AXIS) {
        this.axes[GamepadSticks.AXIS[ax]] = {
            code: GamepadSticks.AXIS[ax],
            name: ax,
            value: 0
        };
    }
};

inherit(GamepadSticks, Input, {
    /**
     * Polls the gamepad object for status updates and emits events when they occur
     *
     * @method pollStatus
     * @param pad {Gamepad} The gamepad object to check
     */
    pollStatus: function(pad) {
        for(var a = 0, al = pad.axes.length; a < al; ++a) {
            var ax = pad.axes[a],
                status = this.axes[a];

            //if we have moved off center by threshold, update the value
            if(Math.abs(ax) >= this.threshold) {
                status.value = ax;
            }
            //otherwise, set it back to zero
            else {
                status.value = 0;
            }

            this.emit(a, status);
        }
    }
});

/**
 * Bindable Gamepad Axes
 *
 * @property AXIS
 * @type Object
 * @static
 */
GamepadSticks.AXIS = {
    LEFT_ANALOGUE_HOR: 0,
    LEFT_ANALOGUE_VERT: 1,
    RIGHT_ANALOGUE_HOR: 2,
    RIGHT_ANALOGUE_VERT: 3
};

GamepadSticks.getGpAxisName = function(i) {
    for(var k in GamepadSticks.AXIS) {
        if(GamepadSticks.AXIS[k] === i) {
            return k;
        }
    }

    return '';
};

module.exports = GamepadSticks;

},{"../../utils/inherit":69,"../Input":39}],45:[function(_dereq_,module,exports){
var Input = _dereq_('../Input'),
    Vector = _dereq_('../../math/Vector'),
    Clock = _dereq_('../../utils/Clock'),
    inherit = _dereq_('../../utils/inherit');

/**
 * Represents a single pointer input method
 *
 * @class Pointer
 * @extends Input
 * @constructor
 * @param id {String|Number} The identifier for this pointer
 * @param manager {Pointers} The pointer manager for this pointer instance
 */
var Pointer = function(id, manager) {
    Input.call(this);

    /**
     * The id of this pointer
     *
     * @property id
     * @type String|Number
     * @readOnly
     */
    this.id = id;

    /**
     * The pointer's manager
     *
     * @property manager
     * @type Pointers
     * @readOnly
     */
    this.manager = manager;

    /**
     * The game instance of the pointer
     *
     * @property game
     * @type Game
     * @readOnly
     */
    this.game = manager.game;

    /**
     * Is this an active pointer (currently touching)?
     *
     * @property active
     * @type Boolean
     * @readOnly
     */
    this.active = false;

    /**
     * Is this the mouse pointer?
     *
     * @property mouse
     * @type Boolean
     * @readOnly
     */
    this.mouse = (id === 1);

    /**
     * The clock for timing stuffz
     *
     * @property clock
     * @type Clock
     * @readOnly
     * @private
     */
    this.clock = new Clock();

    /**
     * The button on the pointer being pressed
     *
     * @property button
     * @type Number
     * @readOnly
     */
    this.button = null;

    /**
     * The type of the pointer
     *
     * @property type
     * @type TYPE
     */
    this.type = null;

    /**
     * Have we emitted the hold event already?
     *
     * @property _holdSent
     * @type Boolean
     * @readOnly
     * @private
     */
    this._holdSent = false;

    this.hitSprites = [];

    this.position = new Vector();
    this.positionDown = new Vector();
};

inherit(Pointer, Input, {
    /**
     * Callback for when a pointerdown event occurs
     *
     * @method down
     * @param evt {DOMEvent} The original DOM Event
     */
    down: function(evt) {
        this.originalEvent = evt;

        //store down timing
        this.timeDown = this.clock.now();
        this.timeHold = 0;
        this._holdSent = false;

        //update x/y position
        this.move(evt);
        this.positionDown.copy(this.position);

        //copy some event vars
        this.button = evt.button;
        this.type = evt.pointerType;

        if(!this.active) {
            this.active = true;
            this.manager.activePointers++;
        }

        for (var i = 0; i < this.hitSprites.length; ++i) {
            var sprite = this.hitSprites[i];

            sprite.__isDown = true;
            sprite.emit('pointerdown', this);
        }
    },
    /**
     * Callback for when a pointerup event occurs
     *
     * @method up
     * @param evt {DOMEvent} The original DOM Event
     */
    up: function(evt) {
        var emit;

        this.originalEvent = evt;

        this.timeUp = this.clock.now();
        this.timeHold = this.timeUp - this.timeDown;

        //consider this a click?
        if(this.timeHold >= 0 && this.timeHold <= this.manager.clickDelay) {
            //is this a double click?
            if((this.timeUp - this.previousClickTime) <= this.manager.doubleClickDelay) {
                emit = 'pointerdoubleclick';
            }
            //only a single click
            else {
                emit = 'pointerclick';
            }

            this.previousClickTime = this.timeUp;
        }

        //mouse is always active
        if(!this.mouse) {
            this.active = false;
            this.manager.activePointers--;
        }

        for(var i = 0; i < this.hitSprites.length; ++i) {
            var sprite = this.hitSprites[i];

            if(sprite.__isDown) {
                sprite.emit('pointerup', this);
                sprite.emit(emit, this);
            }

            sprite.__isDown = false;
        }
    },
    /**
     * Callback for when a pointermove event occurs
     *
     * @method move
     * @param evt {DOMEvent} The original DOM Event
     */
    move: function(evt) {
        this.originalEvent = evt;

        //copy some event vars
        this.button = evt.button;
        this.type = evt.pointerType;

        var rect = this.manager.game.canvas.getBoundingClientRect(); //can we cache this? Maybe update on resize?

        this.position.set(
            (evt.clientX - rect.left) * (this.manager.game.width / rect.width),
            (evt.clientY - rect.top) * (this.manager.game.height / rect.height)
        );

        for(var i = 0; i < this.hitSprites.length; ++i) {
            this.hitSprites[i].emit('pointermove', this);
        }
    },
    /**
     * Callback for when a pointerleave event occurs
     *
     * @method leave
     * @param evt {DOMEvent} The original DOM Event
     */
    leave: function(evt) {
        this.move(evt);

        for(var i = 0; i < this.hitSprites.length; ++i) {
            var sprite = this.hitSprites[i];

            if(sprite.__isOver) {
                sprite.__isOver = false;
                sprite.emit('pointerout', this);
            }
        }
    },
    /**
     * Called internally every frame. Updates the pointer
     *
     * @method update
     * @param dt {Number} The delta time (in seconds) since the last update
     * @private
     */
    update: function() {
        if(!this.active)
            return;

        this.timeHold += this.clock.now() - this.timeDown;

        var holding = (this.timeHold >= this.manager.holdDelay),
            self = this;

        this.checkHits(function(sprite) {
            if(holding && sprite.__isDown) {
                sprite.emit('pointerhold', self);
            }
        });
    },
    hitTest: function(spr) {
        if(!spr.worldVisible)
            return false;

        var worldTransform = spr.worldTransform,
            a00 = worldTransform.a, a01 = worldTransform.b, a02 = worldTransform.tx,
            a10 = worldTransform.c, a11 = worldTransform.d, a12 = worldTransform.ty,
            id = 1 / (a00 * a11 + a01 * -a10),
            x = a11 * id * this.position.x + -a01 * id * this.position.y + (a12 * a01 - a02 * a11) * id,
            y = a00 * id * this.position.y + -a10 * id * this.position.x + (-a12 * a00 + a02 * a10) * id;

        if(spr.hitArea && spr.hitArea.contains) {
            return spr.hitArea.contains(x, y);
        } else {
            var width = spr.texture.frame.width,
                height = spr.texture.frame.height,
                x1 = -width * spr.anchor.x,
                y1;

            if(x > x1 && x < x1 + width) {
                y1 = -height * spr.anchor.y;

                return (y > y1 && y < y1 + height);
            }
        }

        return false;
    },
    checkHits: function(fn) {
        this.hitSprites.length = 0;

        var sprite;

        //hit-test each interactive sprite
        for(var i = 0; i < this.manager.interactiveSprites.length; ++i) {
            sprite = this.manager.interactiveSprites[i];

            if(this.hitTest(sprite)) {
                this.hitSprites.push(sprite);

                if(fn) fn(sprite);

                if(!sprite.__isOver) {
                    sprite.emit('pointerover', this);
                    sprite.__isOver = true;
                }
            } else if(sprite.__isOver) {
                sprite.emit('pointerout', this);
                sprite.__isOver = false;
            }
        }
    },
    getLocalPosition: function(spr) {
        var worldTransform = spr.worldTransform,
            a00 = worldTransform.a, a01 = worldTransform.b, a02 = worldTransform.tx,
            a10 = worldTransform.c, a11 = worldTransform.d, a12 = worldTransform.ty,
            id = 1 / (a00 * a11 + a01 * -a10);

        // set the mouse coords...
        return new Vector(
            a11 * id * this.position.x + -a01 * id * this.position.y + (a12 * a01 - a02 * a11) * id,
            a00 * id * this.position.y + -a10 * id * this.position.x + (-a12 * a00 + a02 * a10) * id
        );
    }
});

/**
 * The type of a pointer
 *
 * @property TYPE
 * @type Object
 */
Pointer.TYPE = {
    TOUCH: 'touch',
    PEN: 'pen',
    MOUSE: 'mouse'
};

module.exports = Pointer;

},{"../../math/Vector":47,"../../utils/Clock":62,"../../utils/inherit":69,"../Input":39}],46:[function(_dereq_,module,exports){
// Thanks to PhotonStorm (http://photonstorm.com/) for this loader!
// heavily insprite by (stolen from): https://github.com/photonstorm/phaser/blob/master/src/loader/Loader.js

var utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    support = _dereq_('../utils/support'),
    EventEmitter = _dereq_('../utils/EventEmitter'),
    C = _dereq_('../constants');

/**
 * The Loader loads and parses different game assets, such as sounds, textures,
 * TMX World files (exported from the [Tiled Editor](http://mapeditor.org)),
 * and Sprite Atlas files (published from [Texture Packer](http://www.codeandweb.com/texturepacker)).
 *
 * @class Loader
 * @extends Object
 * @uses EventEmitter
 * @constructor
 * @param game {Game} Game instance this belongs to
 */
var Loader = function(game) {
    EventEmitter.call(this);

    /**
     * The game instance this loader belongs to
     *
     * @property game
     * @type Game
     */
    this.game = game;

    /**
     * The array of asset keys
     *
     * @property assets
     * @type Array
     */
    this.keys = [];

    /**
     * The asset data
     *
     * @property assets
     * @type Array
     */
    this.assets = {};

    /**
     * Number of assets total to load
     *
     * @property total
     * @type Number
     */
    this.total = 0;

    /**
     * Number of assets done to load (for progress)
     *
     * @property done
     * @type Number
     */
    this.done = 0;

    /**
     * Whether the loader is actively loading the assets
     *
     * @property isLoading
     * @type Boolean
     */
    this.isLoading = false;

    /**
     * Whether the loader has finished loading
     *
     * @property isLoading
     * @type Boolean
     */
    this.hasLoaded = false;

    /**
     * The progress of the loader (0 - 100)
     *
     * @property progress
     * @type Number
     */
    this.progress = 0;

    /**
     * The cross origin value for loading images
     *
     * @property crossOrigin
     * @type String
     */
    this.crossOrigin = '';

    /**
     * The base URL to prepend to a url, requires the trailing slash
     *
     * @property baseUrl
     * @type String
     */
    this.baseUrl = '';

    /**
     * Fired when an item has started loading
     *
     * @event start
     * @param numAssets {Number} The number of assets that are going to be loaded
     */

    /**
     * Fired if a loader encounters an error
     *
     * @event error
     * @param error {mixed} The error that occured when loading
     * @param key {String} The key for the asset that was being loaded
     */

    /**
     * Fired when an item has loaded
     *
     * @event progress
     * @param progress {Number} The integer progress value, between 0 and 100.
     */

    /**
     * Fired when all the assets have loaded
     *
     * @event complete
     */
};

inherit(Loader, Object, {
    /**
     * Check whether asset exists with a specific key.
     *
     * @method hasKey
     * @param key {String} Key of the asset you want to check.
     * @return {Boolean} Return true if exists, otherwise return false.
     */
    hasKey: function(key) {
        return !!this.assets[key];
    },

    /**
     * Reset loader, this will remove all loaded assets from the loader's stored list (but not from the cache).
     *
     * @method reset
     * @return {Loader} Returns itself.
     * @chainable
     */
    reset: function() {
        this.progress = 0;
        this.total = 0;
        this.done = 0;
        this.hasLoaded = false;
        this.isLoading = false;
        this.assets = {};
        this.keys.length = 0;

        return this;
    },

    /**
     * Adds an asset to be loaded
     *
     * @method add
     * @param type {String} The type of asset ot load (image, spritesheet, textureatlas, bitmapfont, tilemap, tileset, audio, etc)
     * @param key {String} The unique key of the asset to identify it
     * @param url {String} The URL to load the resource from
     * @param [options] {Object} Extra options to apply to the asset, different asset types may require extra options
     * @param [options.crossOrigin=false] {Boolean} True if an image load should be treated as crossOrigin
     * @return {Loader} Returns itself.
     * @chainable
     */
    add: function(type, key, url, opts) {
        var entry = {
            type: type,
            key: key,
            url: url,
            image: null,
            data: null,
            error: false,
            loaded: false
        };

        if(opts !== undefined) {
            for(var p in opts) {
                entry[p] = opts[p];
            }
        }

        this.assets[key] = entry;
        this.keys.push(key);
        this.total++;

        return this;
    },

    /**
     * Add an image to the Loader.
     *
     * @method image
     * @param key {String} Unique asset key of this image file.
     * @param url {String} URL of image file.
     * @param [overwrite=false] {Boolean} If an entry with a matching key already exists this will over-write it.
     * @return {Loader} Returns itself.
     * @chainable
     */
    image: function(key, url, overwrite) {
        if(overwrite || !this.hasKey(key))
            this.add('image', key, url);

        return this;
    },

    /**
     * Add a text file to the Loader.
     *
     * @method text
     * @param key {String} Unique asset key of this image file.
     * @param url {String} URL of image file.
     * @param [overwrite=false] {Boolean} If an entry with a matching key already exists this will over-write it.
     * @return {Loader} Returns itself.
     * @chainable
     */
    text: function(key, url, overwrite) {
        if(overwrite || !this.hasKey(key))
            this.add('text', key, url);

        return this;
    },

    /**
     * Add a sprite sheet image to the Loader.
     *
     * @method spritesheet
     * @param key {String} Unique asset key of this image file.
     * @param url {String} URL of image file.
     * @param frameWidth {Number} Width of each single frame.
     * @param frameHeight {Number} Height of each single frame.
     * @param numFrames {Number} How many frames in this sprite sheet.
     * @param [overwrite=false] {Boolean} If an entry with a matching key already exists this will over-write it.
     * @return {Loader} Returns itself.
     * @chainable
     */
    spritesheet: function(key, url, frameWidth, frameHeight, numFrames, overwrite) {
        if(overwrite || !this.hasKey(key))
            this.add('spritesheet', key, url, {
                frameWidth: frameWidth,
                frameHeight: frameHeight,
                numFrames: numFrames
            });

        return this;
    },

    /**
     * Add an audio file to the Loader.
     *
     * @method audio
     * @param key {String} Unique asset key of this image file.
     * @param urls {Array<String>} URLs of audio files.
     * @param [overwrite=false] {Boolean} If an entry with a matching key already exists this will over-write it.
     * @return {Loader} Returns itself.
     * @chainable
     */
    audio: function(key, urls, overwrite) {
        if(overwrite || !this.hasKey(key))
            this.add('audio', key, urls);

        return this;
    },

    /**
     * Add a tilemap to the Loader.
     *
     * @method tilemap
     * @param key {String} Unique asset key of the tilemap data.
     * @param url {String} The url of the map data file (csv/json/xml)
     * @param [data] {String|Object} The data for the map, (to use instead of loading from a URL)
     * @param [format=FILE_FORMAT.JSON] {Number} The format of the map data.
     * @param [overwrite=false] {Boolean} If an entry with a matching key already exists this will over-write it.
     * @return {Loader} Returns itself.
     * @chainable
     */
    tilemap: function(key, url, data, format, overwrite) {
        if(overwrite || !this.hasKey(key)) {
            if(!format) format = C.FILE_FORMAT.JSON;

            if(typeof data === 'string') {
                switch(format) {
                    case C.FILE_FORMAT.JSON:
                        data = JSON.parse(data);
                        break;

                    case C.FILE_FORMAT.XML:
                        data = C.utils.parseXML(data);
                        break;

                    case C.FILE_FORMAT.CSV:
                        break;
                }
            }

            this.add('tilemap', key, url, {
                data: data,
                format: format
            });
        }

        return this;
    },

    /**
     * Add a bitmap font to the Loader.
     *
     * @method bitmapFont
     * @param key {String} Unique asset key of the bitmap font.
     * @param textureURL {String} The url of the font image file.
     * @param [dataUrl] {String} The url of the font data file (xml/fnt)
     * @param [data] {Object} An optional XML data object (to use instead of loading from a URL)
     * @param [format=FILE_FORMAT.XML] {FILE_FORMAT} The format of the bitmap font data.
     * @param [overwrite=false] {Boolean} If an entry with a matching key already exists this will over-write it.
     * @return {Loader} Returns itself.
     * @chainable
     */
    bitmapFont: function(key, textureUrl, dataUrl, data, format, overwrite) {
        if(overwrite || !this.hasKey(key)) {
            if(!format) format = C.FILE_FORMAT.XML;

            if(typeof data === 'string') {
                switch(format) {
                    case C.FILE_FORMAT.XML:
                        data = utils.parseXML(data);
                        break;

                    case C.FILE_FORMAT.JSON:
                        data = JSON.parse(data);
                        break;
                }
            }

            this.add('bitmapfont', key, textureUrl, {
                dataUrl: dataUrl,
                data: data,
                format: format
            });
        }

        return this;
    },

    /**
     * Add a JSON-Array formatted texture atlas. Equivalent to running
     * `atlas(key, textureURL, dataUrl, data, gf.ATLAS_FORMAT.JSON_ARRAY);`
     *
     * @param key {string} Unique asset key of the texture atlas file.
     * @param textureUrl {string} The url of the texture atlas image file.
     * @param [dataUrl] {string} The url of the texture atlas data file (json/xml)
     * @param [data] {object} A JSON or XML data object (to use instead of loading from a URL)
     * @return {Loader} Returns itself.
     * @chainable
     */
    atlasJSONArray: function(key, textureURL, dataUrl, data) {
        return this.atlas(key, textureURL, dataUrl, data, C.ATLAS_FORMAT.JSON_ARRAY);
    },

    /**
     * Add a JSON-Hash formatted texture atlas. Equivalent to running
     * `atlas(key, textureURL, dataUrl, data, gf.ATLAS_FORMAT.JSON_HASH);`
     *
     * @param key {string} Unique asset key of the texture atlas file.
     * @param textureUrl {string} The url of the texture atlas image file.
     * @param [dataUrl] {string} The url of the texture atlas data file (json/xml)
     * @param [data] {object} A JSON or XML data object (to use instead of loading from a URL)
     * @return {Loader} Returns itself.
     * @chainable
     */
    atlasJSONHash: function(key, textureURL, dataUrl, data) {
        return this.atlas(key, textureURL, dataUrl, data, C.ATLAS_FORMAT.JSON_HASH);
    },

    /**
     * Add an XML formatted texture atlas. Equivalent to running
     * `atlas(key, textureURL, dataUrl, data, gf.ATLAS_FORMAT.XML_STARLING);`
     *
     * @param key {string} Unique asset key of the texture atlas file.
     * @param textureUrl {string} The url of the texture atlas image file.
     * @param [dataUrl] {string} The url of the texture atlas data file (json/xml)
     * @param [data] {object} A JSON or XML data object (to use instead of loading from a URL)
     * @return {Loader} Returns itself.
     * @chainable
     */
    atlasXML: function(key, textureURL, dataUrl, data) {
        return this.atlas(key, textureURL, dataUrl, data, C.ATLAS_FORMAT.XML_STARLING);
    },

    /**
     * Add a new texture atlas loading request.
     * @param key {string} Unique asset key of the texture atlas file.
     * @param textureUrl {string} The url of the texture atlas image file.
     * @param [dataUrl] {string} The url of the texture atlas data file (json/xml)
     * @param [data] {object} A JSON or XML data object (to use instead of loading from a URL)
     * @param [format] {number} A value describing the format of the data.
     * @param [overwrite=false] {Boolean} If an entry with a matching key already exists this will over-write it.
     * @return {Loader} Returns itself.
     * @chainable
     */
    atlas: function(key, textureUrl, dataUrl, data, format, overwrite) {
        if(overwrite || !this.hasKey(key)) {
            if(!format) format = C.ATLAS_FORMAT.JSON_ARRAY;

            if(typeof data === 'string') {
                switch(format) {
                    case C.ATLAS_FORMAT.XML_STARLING:
                        data = utils.parseXML(data);
                        break;

                    case C.ATLAS_FORMAT.JSON_ARRAY:
                    case C.ATLAS_FORMAT.JSON_HASH:
                        data = JSON.parse(data);
                        break;
                }
            }

            this.add('textureatlas', key, textureUrl, {
                dataUrl: dataUrl,
                data: data,
                format: format
            });
        }

        return this;
    },

    /**
     * Starts the loading of all the assets that are queued to load
     *
     * @method start
     * @return {Loader} Returns itself.
     * @chainable
     */
    start: function() {
        if(this.isLoading) return;

        this.progress = 0;
        this.hasLoaded = false;
        this.isLoading = true;

        this.emit('start', this.keys.length);

        if(this.keys.length > 0) {
            while(this.keys.length > 0)
                this.loadFile();
        } else {
            this.isLoading = false;
            this.progress = 100;
            this.hasLoaded = true;
            this.emit('complete');
        }

        return this;
    },

    /**
     * Loads a single asset from the queued assets in this Loader. To load a single file first queue it by using
     * one of the methods named for an asset (like `audio`, `image`, `tilemap`, etc.), then call this to load the
     * first in the queue.
     *
     * Note: To load the entire queue at once use `start`.
     *
     * @method loadFile
     * @return {Loader} Returns itself.
     * @chainable
     */
    loadFile: function() {
        var file = this.assets[this.keys.shift()],
            self = this;

        switch(file.type) {
            //load images
            case 'image':
            case 'spritesheet':
            case 'textureatlas':
            case 'bitmapfont':
                file.image = new Image();
                file.image.name = file.key;
                file.image.addEventListener('load', this.fileComplete.bind(this, file.key), false);
                file.image.addEventListener('error', this.fileError.bind(this, file.key), false);
                file.image.crossOrigin = file.crossOrigin !== undefined ? file.crossOrigin : this.crossOrigin;
                file.image.src = this.baseUrl + file.url;
                break;

            //load tilemap
            case 'tilemap':
                utils.ajax({
                    url: this.baseUrl + file.url,
                    dataType: this._getFormatAjaxType(file.format),
                    load: function(data) {
                        file.data = data;
                        self.fileComplete(file.key);
                    },
                    error: function(err) {
                        self.fileError(file.key, err);
                    }
                });
                break;

            //load audio
            case 'audio':
                file.url = this.getAudioUrl(file.url);

                if(file.url) {
                    if(support.webAudio) {
                        utils.ajax({
                            url: this.baseUrl + file.url,
                            dataType: 'arraybuffer',
                            load: function(data) {
                                file.data = data;
                                self.fileComplete(file.key);
                            },
                            error: function(err) {
                                self.fileError(file.key, err);
                            }
                        });
                    } else if(support.htmlAudio) {
                        file.data = new Audio();
                        file.data.name = file.key;
                        file.data.preload = 'auto';
                        file.data.src = this.baseUrl + file.url;

                        file.data.addEventListener('error', file._bndError = this.fileError.bind(this, file.key), false);
                        file.data.addEventListener('canplaythrough', file._bndComplete = this.fileComplete.bind(this, file.key), false);
                        file.data.load();
                    }
                } else {
                    this.fileError(file.key, 'No supported audio URL could be determined!');
                }
                break;

            case 'text':
                utils.ajax({
                    url: this.baseUrl + file.url,
                    dataType: 'text',
                    load: function(data) {
                        file.data = data;
                        self.fileComplete(file.key);
                    },
                    error: function(err) {
                        self.fileError(file.key, err);
                    }
                });
                break;
        }

        return this;
    },

    /**
     * Chooses the audio url to use based on browser support.
     *
     * @method getAudioUrl
     * @param urls {Array<String>} An array of URLs to choose from, chooses the first in the array to be
     *      supported by the browser.
     * @return {String} Returns the URL that was chosen, or `undefined` if none are supported.
     */
    getAudioUrl: function(urls) {
        for(var i = 0, il = urls.length; i < il; ++i) {
            var url = urls[i],
                ext = url.match(/.+\.([^?]+)(\?|$)/);

            ext = (ext && ext.length >= 2) ? ext[1] : url.match(/data\:audio\/([^?]+);/)[1];

            //if we can play this url, then set the source of the player
            if(support.codec && support.codec[ext]) {
                return url;
            }
        }

        return false;
    },

    /**
     * Error occured when load a file.
     *
     * @method fileError
     * @param key {String} Key of the error loading file.
     * @param error {mixed} The error that was thrown.
     * @private
     */
    fileError: function(key, error) {
        this.assets[key].loaded = true;
        this.assets[key].error = error;

        this.fileDone(key, error);
    },

    /**
     * Called when a file is successfully loaded.
     *
     * @method fileComplete
     * @param key {string} Key of the successfully loaded file.
     * @private
     */
    fileComplete: function(key) {
        if(!this.assets[key])
            return utils.warn('fileComplete key is invalid!', key);

        this.assets[key].loaded = true;

        var file = this.assets[key],
            done = true,
            self = this;

        switch(file.type) {
            case 'image':
                this.game.cache.addImage(file);
                break;

            case 'spritesheet':
                this.game.cache.addSpriteSheet(file);
                break;

            case 'tilemap':
                file.baseUrl = file.url.replace(/[^\/]*$/, '');
                file.numImages = file.numLoaded = 0;
                file.images = {};

                if(file.format === C.FILE_FORMAT.JSON) {
                    done = false;
                    this._loadJsonTilesets(file);
                } else if(file.format === C.FILE_FORMAT.XML) {
                    done = false;
                    this._loadXmlTilesets(file);
                }
                break;

            case 'textureatlas':
                done = false;
                this._dataget(file, function() {
                    self.game.cache.addTextureAtlas(file);
                });
                break;

            case 'bitmapfont':
                done = false;
                this._dataget(file, function() {
                    self.game.cache.addBitmapFont(file);
                });
                break;

            case 'audio':
                if(support.webAudio) {
                    file.webAudio = true;
                    file.decoded = false;
                } else {
                    file.data.removeEventListener('error', file._bndError);
                    file.data.removeEventListener('canplaythrough', file._bndComplete);
                }

                this.game.cache.addAudio(file);
                break;

            case 'text':
                this.game.cache.addText(file);
                break;
        }

        if(done) {
            this.fileDone(file.key);
        }
    },

    /**
     * Called when a file is done (error or loaded)
     *
     * @method fileDone
     * @param key {String} Key of the file done
     * @param error {mixed} The error that occurred (if there was one)
     * @private
     */
    fileDone: function(key, error) {
        this.done++;
        this.progress = Math.round((this.done / this.total) * 100);

        this.emit('progress', this.progress);

        if(error) {
            utils.warn('Error loading file "' + key + '", error received:', error);
            this.emit('error', error, key);
        }

        if(this.progress >= 100) {
            this.progress = 100;
            this.hasLoaded = true;
            this.isLoading = false;

            this.emit('complete');
        }
    },

    /**
     * Returns the ajax type that represents each format type
     *
     * @method _getFormatAjaxType
     * @param type {ATLAS_FORMAT|FILE_FORMAT} The format to get an ajax type for
     * @private
     */
    _getFormatAjaxType: function(type) {
        switch(type) {
            case C.ATLAS_FORMAT.JSON_ARRAY:
            case C.ATLAS_FORMAT.JSON_HASH:
            case C.FILE_FORMAT.JSON:
                return 'json';

            case C.ATLAS_FORMAT.XML_STARLING:
            case C.FILE_FORMAT.XML:
                return 'xml';

            case C.FILE_FORMAT.CSV:
                return 'text';
        }
    },

    /**
     * Gets a file's data via ajax.
     *
     * @method _dataget
     * @param file {Object} The file descriptor object
     * @param [callback] {Function} The callback to call once the file has loaded. `fileDone` or `fileError` will be
     *      called for you.
     * @private
     */
    _dataget: function(file, cb) {
        var self = this;

        if(!file.dataUrl) {
            setTimeout(cb, 1);
        } else {
            utils.ajax({
                url: this.baseUrl + file.dataUrl,
                dataType: this._getFormatAjaxType(file.format),
                load: function(data) {
                    file.data = data;
                    if(cb) cb();
                    self.fileDone(file.key);
                },
                error: function(err) {
                    self.fileError(file.key, err);
                }
            });
        }
    },

    /**
     * Loads the tilesets found in a JSON formatted tilemap object.
     *
     * @method _loadJsonTilesets
     * @param file {Object} The file descriptor object
     * @private
     */
    _loadJsonTilesets: function(file) {
        var data = file.data;

        //loop through each tileset and load the texture
        for(var i = 0, il = data.tilesets.length; i < il; ++i) {
            var set = data.tilesets[i];

            if(set.image) {
                file.numImages++;
                file.images[set.image] = this._getTilesetImage(file, set.image);
            } else if(set.tiles) {
                for(var t in set.tiles) {
                    file.numImages++;
                    file.images[set.tiles[t].image] = this._getTilesetImage(file, set.tiles[t].image);
                }
            }
        }
    },

    /**
     * Loads the tilesets found in a XML formatted tilemap object.
     *
     * @method _loadXmlTilesets
     * @param file {Object} The file descriptor object
     * @private
     */
    _loadXmlTilesets: function(file) {
        var data = file.data,
            tilesets = data.getElementsByTagName('tileset');

        for(var i = 0, il = tilesets.length; i < il; ++i) {
            var set = tilesets[i],
                imgElm = set.getElementsByTagName('image')[0];

            if(!imgElm) continue;

            file.numImages++;

            var src = imgElm.attributes.getNamedItem('source').nodeValue;
            file.images[src] = this._getTilesetImage(file, src);
        }
    },

    _getTilesetImage: function(file, src) {
        var img = new Image();

        img.addEventListener('load', this._onTilesetLoaded.bind(this, file), false);
        img.addEventListener('error', this._onTilesetError.bind(this, file), false);
        img.crossOrigin = file.crossOrigin !== undefined ? file.crossOrigin : this.crossOrigin;
        img.src = this.baseUrl + file.baseUrl + src;

        return img;
    },

    /**
     * Called each time a tileset is loaded successfully.
     *
     * @method _onTilesetLoaded
     * @param file {Object} The file descriptor object.
     * @private
     */
    _onTilesetLoaded: function(file) {
        file.numLoaded++;

        if(file.numImages === file.numLoaded) {
            this.game.cache.addTilemap(file);
            this.fileDone(file.key);
        }
    },

    /**
     * Called each time a tileset has an error when loading.
     *
     * @method _onTilesetError
     * @param file {Object} The file descriptor object.
     * @param error {mixed} The error thrown when loading.
     * @private
     */
    _onTilesetError: function(file, error) {
        file.error = error;
        file.numLoaded++;

        if(file.numImages === file.numLoaded) {
            this.fileDone(file.key, error);
        }
    }
});

module.exports = Loader;

},{"../constants":11,"../utils/EventEmitter":64,"../utils/inherit":69,"../utils/support":70,"../utils/utils":71}],47:[function(_dereq_,module,exports){
var inherit = _dereq_('../utils/inherit'),
    math = _dereq_('./math');

/**
 * A 2d Vector implementation stolen directly from mrdoob's THREE.js
 * [Vector2d](https://github.com/mrdoob/three.js/blob/master/src/math/Vector2.js)
 *
 * @class Vector
 * @extends Object
 * @constructor
 * @param x {Number} The x component of the vector
 * @param y {Number} The y component of the vector
 */
var Vector = function(x, y) {
    this.x = x || 0;
    this.y = y || 0;
};

inherit(Vector, Object, {
    /**
     * Sets the value of the vector
     *
     * @method set
     * @param x {Number} The x component of the vector
     * @param y {Number} The y component of the vector
     * @return {Vector} Returns itself.
     * @chainable
     */
    set: function(x, y) {
        this.x = x;
        this.y = y;

        return this;
    },
    /**
     * Sets the X value of the vector
     *
     * @method setX
     * @param x {Number} The x component of the vector
     * @return {Vector} Returns itself.
     * @chainable
     */
    setX: function(x) {
        this.x = x;

        return this;
    },
    /**
     * Sets the Y value of the vector
     *
     * @method setY
     * @param y {Number} The y component of the vector
     * @return {Vector} Returns itself.
     * @chainable
     */
    setY: function(y) {
        this.y = y;

        return this;
    },
    /**
     * Sets a component value of the vector
     *
     * @method setComponent
     * @param index {Number} The index of the component to set (0 = x, 1 = y)
     * @param value {Number} The value to set the component to
     * @return {Vector} Returns itself.
     * @chainable
     */
    setComponent: function(index, value) {
        switch(index) {
            case 0:
                this.x = value;
                break;
            case 1:
                this.y = value;
                break;
            default:
                throw new RangeError('index is out of range: ' + index);
        }

        return this;
    },
    /**
     * Gets a component value of the vector
     *
     * @method getComponent
     * @param index {Number} The index of the component to set (0 = x, 1 = y)
     * @return {Number} Returns the component value
     */
    getComponent: function(index) {
        switch(index) {
            case 0:
                return this.x;
            case 1:
                return this.y;
            default:
                throw new RangeError('index is out of range: ' + index);
        }
    },
    /**
     * Copies the passed vector's components to this vector
     *
     * @method copy
     * @param vector {Vector} The vector to copy the values from
     * @return {Vector} Returns itself.
     * @chainable
     */
    copy: function(v) {
        this.x = v.x;
        this.y = v.y;

        return this;
    },
    /**
     * Floors the vector components
     *
     * @method floor
     * @return {Vector} Returns itself.
     * @chainable
     */
    floor: function () {
        this.x = math.floor(this.x);
        this.y = math.floor(this.y);

        return this;
    },
    /**
     * Ceils the vector components
     *
     * @method ceil
     * @return {Vector} Returns itself.
     * @chainable
     */
    ceil: function () {
        this.x = math.ceil(this.x);
        this.y = math.ceil(this.y);

        return this;
    },
    round: function() {
        this.x = math.round(this.x);
        this.y = math.round(this.y);

        return this;
    },
    /**
     * Adds a vector to this one
     *
     * @method add
     * @param vector {Vector} The vector to add to this one
     * @return {Vector} Returns itself.
     * @chainable
     */
    add: function(v) {
        this.x += v.x;
        this.y += v.y;

        return this;
    },
    /**
     * Adds two vectors to each other and stores the result in this vector
     *
     * @method addVectors
     * @param vector1 {Vector}
     * @param vector2 {Vector}
     * @return {Vector} Returns itself.
     * @chainable
     */
    addVectors: function(a, b) {
        this.x = a.x + b.x;
        this.y = a.y + b.y;

        return this;
    },
    /**
     * Adds a scalar value to the x and y components of this vector
     *
     * @method addScalar
     * @param scalar {Number} The scalar value to add
     * @return {Vector} Returns itself.
     * @chainable
     */
    addScalar: function(s) {
        this.x += s;
        this.y += s;

        return this;
    },
    /**
     * Subtracts a vector from this one
     *
     * @method sub
     * @param vector {Vector} The vector to subtract from this one
     * @return {Vector} Returns itself.
     * @chainable
     */
    sub: function(v) {
        this.x -= v.x;
        this.y -= v.y;

        return this;
    },
    /**
     * Subtracts two vectors from each other and stores the result in this vector
     *
     * @method subVectors
     * @param vector1 {Vector}
     * @param vector2 {Vector}
     * @return {Vector} Returns itself.
     * @chainable
     */
    subVectors: function(a, b) {
        this.x = a.x - b.x;
        this.y = a.y - b.y;

        return this;
    },
    /**
     * Multiplies the x and y components of this vector by a scalar value
     *
     * @method multiplyScalar
     * @param scalar {Number} The value to multiply by
     * @return {Vector} Returns itself.
     * @chainable
     */
    multiplyScalar: function(s) {
        this.x *= s;
        this.y *= s;

        return this;
    },
    /**
     * Divides the x and y components of this vector by a scalar value
     *
     * @method divideScalar
     * @param scalar {Number} The value to divide by
     * @return {Vector} Returns itself.
     * @chainable
     */
    divideScalar: function(s) {
        if(s !== 0) {
            this.x /= s;
            this.y /= s;
        } else {
            this.set(0, 0);
        }

        return this;
    },
    /**
     * Sets this vector components to the minimum value when compared to the passed vector's components
     *
     * @method min
     * @param vector {Vector} The vector to compare to
     * @return {Vector} Returns itself.
     * @chainable
     */
    min: function(v) {
        if(this.x > v.x) {
            this.x = v.x;
        }

        if(this.y > v.y) {
            this.y = v.y;
        }

        return this;
    },
    /**
     * Sets this vector components to the maximum value when compared to the passed vector's components
     *
     * @method max
     * @param vector {Vector} The vector to compare to
     * @return {Vector} Returns itself.
     * @chainable
     */
    max: function(v) {
        if(this.x < v.x) {
            this.x = v.x;
        }

        if(this.y < v.y) {
            this.y = v.y;
        }

        return this;
    },
    /**
     * Clamps the vectors components to be between min and max
     *
     * @method max
     * @param min {Vector} The minimum value a component can be
     * @param max {Vector} The maximum value a component can be
     * @return {Vector} Returns itself.
     * @chainable
     */
    clamp: function(min, max) {
        // This function assumes min < max, if this assumption
        //isn't true it will not operate correctly
        if(this.x < min.x) {
            this.x = min.x;
        } else if(this.x > max.x) {
            this.x = max.x;
        }

        if(this.y < min.y) {
            this.y = min.y;
        } else if(this.y > max.y) {
            this.y = max.y;
        }

        return this;
    },
    /**
     * Negates this vector (multiplies by -1)
     *
     * @method negate
     * @return {Vector} Returns itself.
     * @chainable
     */
    negate: function() {
        return this.multiplyScalar(-1);
    },
    /**
     * Project this vector on to another vector.
     *
     * @param v {Vector} The vector to project onto.
     * @return {Vector} Returns itself.
     * @chainable
     */
    project: function(v) {
        var amt = this.dot(v) / v.lengthSq();
        this.x = amt * v.x;
        this.y = amt * v.y;

        return this;
    },
    /**
     * Project this vector onto a vector of unit length.
     *
     * @param v {Vector} The unit vector to project onto.
     * @return {Vector} Returns itself.
     * @chainable
     */
    projectN: function(v) {
        var amt = this.dot(v);
        this.x = amt * v.x;
        this.y = amt * v.y;

        return this;
    },
    /**
     * Reflect this vector on an arbitrary axis.
     *
     * @param axis {Vector} The vector representing the axis.
     * @return {Vector} Returns itself.
     * @chainable
     */
    reflect: function(axis) {
        var x = this.x;
        var y = this.y;
        this.project(axis).multiplyScalar(2);
        this.x -= x;
        this.y -= y;

        return this;
    },
    /**
     * Reflect this vector on an arbitrary axis (represented by a unit vector)
     *
     * @param axis {Vector} The unit vector representing the axis.
     * @return {Vector} Returns itself.
     * @chainable
     */
    reflectN: function(axis) {
        var x = this.x;
        var y = this.y;
        this.projectN(axis).multiplyScalar(2);
        this.x -= x;
        this.y -= y;

        return this;
    },
    /**
     * Performs the dot product between this vector and the passed one and returns the result
     *
     * @method dot
     * @param vector {Vector}
     * @return {Number} Returns the dot product
     */
    dot: function(v) {
        return this.x * v.x + this.y * v.y;
    },
    /**
     * Calculates the square length of the vector
     *
     * @method lengthSq
     * @return {Number} Returns the square length of the vector
     */
    lengthSq: function() {
        return this.dot(this);
    },
    /**
     * Calculates the length of the vector
     *
     * @method length
     * @return {Number} Returns the length of the vector
     */
    length: function() {
        return math.sqrt(this.lengthSq());
    },
    /**
     * Normalizes this vector (divides by its length)
     *
     * @method normalize
     * @return {Vector} Returns the normalized vector
     */
    normalize: function() {
        return this.divideScalar(this.length());
    },
    /**
     * Calculates the distance to the passed vector
     *
     * @method distanceTo
     * @param vector {Vector} The vector to check distance to
     * @return {Number} The distance
     */
    distanceTo: function(v) {
        return math.sqrt(this.distanceToSquared(v));
    },
    /**
     * Calculates the square distance to the passed vector
     *
     * @method distanceToSquared
     * @param vector {Vector} The vector to check distance to
     * @return {Number} The square distance
     */
    distanceToSquared: function(v) {
        var dx = this.x - v.x, dy = this.y - v.y;
        return dx * dx + dy * dy;
    },
    /**
     * Sets the length of the vector
     *
     * @method setLength
     * @param length {Number} The length to set this vector to
     * @return {Vector} Returns itself.
     * @chainable
     */
    setLength: function(l) {
        var oldLength = this.length();

        if(oldLength !== 0 && l !== oldLength) {
            this.multiplyScalar(l / oldLength);
        }

        return this;
    },
    /**
     * Performs a linear interpolation between this vector and the passed vector
     *
     * @method lerp
     * @param vector {Vector} The vector to interpolate with
     * @param alpha {Number} The amount to interpolate [0-1] or extrapolate (1-]
     * @return {Vector} Returns itself.
     * @chainable
     */
    lerp: function(v, alpha) {
        this.x += (v.x - this.x) * alpha;
        this.y += (v.y - this.y) * alpha;

        return this;
    },
    /**
     * Rotates the vector by 90 degrees
     *
     * @return {Vector} Returns itself.
     * @chainable
     */
    perp: function() {
        var x = this.x;
        this.x = this.y;
        this.y = -x;

        return this;
    },
    /**
     * Rotates the vector by an arbitrary angle around an arbitrary point in space
     *
     * @method rotate
     * @param angle {Number} The angle in radians to rotate by
     * @param anchor {Vector} The anchor point to rotate around
     * @return {Vector} Returns itself.
     * @chainable
     */
    rotate: function(angle, anchor) {
        var dist = anchor.distanceTo(this);

        return this.set(
            anchor.x + (dist * math.cos(angle)),
            anchor.y + (dist * math.sin(angle))
        );
    },
    /**
     * Checks if this vector is equal to another
     *
     * @method equals
     * @param vector {Vector} The vector to compare with
     * @return {Vector} Returns itself.
     * @chainable
     */
    equals: function(v) {
        return ((v.x === this.x) && (v.y === this.y));
    },
    /**
     * Returns an array with the components of this vector as the elements
     *
     * @method toArray
     * @return {Vector} Returns an array of [x,y] form
     */
    toArray: function () {
        return [this.x, this.y];
    },
    /**
     * Creates a new instance of Vector, with the same components as this vector
     *
     * @method clone
     * @return {Vector} Returns a new Vector with the same values
     */
    clone: function () {
        return new Vector(this.x, this.y);
    }
});

/**
 * A vector that is always 0,0
 *
 * @property ZERO
 * @type Vector
 * @readOnly
 * @static
 * @final
 */
Vector.ZERO = new Vector();

module.exports = Vector;

},{"../utils/inherit":69,"./math":48}],48:[function(_dereq_,module,exports){
var random = _dereq_('./random'),
    PIXI = _dereq_('pixi.js');

/**
 * The grapefruit math library, used to abstract commonly used math operations
 *
 * @class math
 * @extends Object
 * @static
 */
var math = {
    /**
     * The factor to multiply by to convert Degrees into Radians. The value is π/180
     *
     * @property DEG_TO_RAD
     * @type Number
     * @readOnly
     */
    DEG_TO_RAD: Math.PI / 180,
    /**
     * The factor to multiply by to convert Radians into Degrees. The value is 180/π
     *
     * @property RAD_TO_DEG
     * @type Number
     * @readOnly
     */
    RAD_TO_DEG: 180 / Math.PI,
    /**
     * Contains all the functions for generating deterministic random values.
     *
     * @property rand
     * @type random
     * @readOnly
     */
    rand: random,
    /**
     * A Matrix class, directory exposes PIXI.Matrix.
     *
     * @property Matrix
     * @type Matrix
     */
    Matrix: PIXI.Matrix,

    /**
     * Alias some native functions for great justice (or incase we want to override)
     */

    /**
     * Applys a Floor operation to a value, currently uses native Math.floor
     * since quicker solutions like `~~value` or `value | 0` only deal with 32-bits.
     * For example `~~760895687099.0011` is `686475707` which is wrong.
     *
     * @method floor
     * @param num {Number} The number to floor
     * @return {Number} The floored value
     */
    floor: Math.floor,
    /**
     * Applys a Ceiling operation to a value, currently uses native Math.ceil
     * since it deals with all edge cases
     *
     * @method ceil
     * @param num {Number} The number to ceil
     * @return {Number} The ceiling value
     */
    ceil: Math.ceil,
    /**
     * Returns the absolute value of a number, currently uses native Math.abs
     * since it is more performant than tricks you can use.
     * see:
     *      http://jsperf.com/math-abs-vs-bitwise/7
     *      http://jsperf.com/abs-value
     *      http://jsperf.com/math-abs-vs-bitwise/3
     *
     * @method abs
     * @param num {Number} The number to get the absolute value for
     * @return {Number} The absolute value
     */
    abs: Math.abs,
    /**
     * Returns the square root of a number, currently uses native Math.sqrt
     *
     * @method sqrt
     * @param num {Number} The number to get the sqrt of
     * @return {Number} The sqrt value
     */
    sqrt: Math.sqrt,
    /**
     * Returns the min of the values passed, currently uses native Math.min
     *
     * @method min
     * @param num* {Number...} The numbers to compare
     * @return {Number} The min value
     */
    min: Math.min,
    /**
     * Returns the max of the values passed, currently uses native Math.max
     *
     * @method max
     * @param num* {Number...} The numbers to compare
     * @return {Number} The max value
     */
    max: Math.max,
    /**
     * Rounds a number to the closest integer value (0.5 goes up), currently
     * uses native Math.round since in modern browsers it is faster that the
     * different tricks and will operate in the proper bit width.
     *
     * @method round
     * @param num {Number} The number to round
     * @return {Number} The rounded value
     */
    round: Math.round,

    /**
     * Clamps a number between two values.
     *
     * @method clamp
     * @param num {Number} The number to clamp
     * @param min {Number} The minimum value the number is allowed to be
     * @param max {Number} The maximum value the number is allowed to be
     * @return {Number} The clamped value
     */
    clamp: function(n, min, max) {
        return math.max(min, math.min(max, n));
    },
    /**
     * Truncates the decimal from a number
     *
     * @method truncate
     * @param num {Number} The number to truncate
     * @return {Number} The truncated value
     */
    truncate: function(n) {
        return (n > 0) ? math.floor(n) : math.ceil(n);
    },
    /**
     * Snaps a number to a grid value.
     * For example, if you have a grid with gaps the size of 10 horizontally, and
     * a position of 11, it would snap to 10; a position of 18 would snap to 20
     *
     * @method snap
     * @param num {Number} The number to snap
     * @param gap {Number} The gap size of the grid (the tile size)
     * @param [offset=0] {Number} The starting offset of a grid slice (aka tile)
     * @return {Number} The snapped value
     */
    snap: function(n, gap, offset) {
        if(gap === 0) return n;

        offset = offset || 0;

        n -= offset;
        n = gap * math.round(n / gap);

        return offset + n;
    },
    /**
     * Snaps a number to a grid value, using floor.
     * For example, if you have a grid with gaps the size of 10 horizontally, and
     * a position of 11, it would snap to 10; a position of 18 would also snap to 10
     *
     * @method snapFloor
     * @param num {Number} The number to snap
     * @param gap {Number} The gap size of the grid (the tile size)
     * @param [offset=0] {Number} The starting offset of a grid slice (aka tile)
     * @return {Number} The snapped value
     */
    snapFloor: function(n, gap, offset) {
        if(gap === 0) return n;

        offset = offset || 0;

        n -= offset;
        n = gap * math.floor(n / gap);

        return offset + n;
    },
    /**
     * Snaps a number to a grid value, using ceiling.
     * For example, if you have a grid with gaps the size of 10 horizontally, and
     * a position of 11, it would snap to 20; a position of 18 would also snap to 20
     *
     * @method snapCeil
     * @param num {Number} The number to snap
     * @param gap {Number} The gap size of the grid (the tile size)
     * @param [offset=0] {Number} The starting offset of a grid slice (aka tile)
     * @return {Number} The snapped value
     */
    snapCeil: function(n, gap, offset) {
        if(gap === 0) return n;

        offset = offset || 0;

        n -= offset;
        n = gap * math.ceil(n / gap);

        return offset + n;
    },
    /**
     * Convert radians to degrees
     *
     * @method radiansToDegrees
     * @param angle {Number} The angle in radians to convert
     * @return {Number} The angle in degrees
     */
    radiansToDegrees: function(angle) {
        return angle * math.RAD_TO_DEG;
    },
    /**
     * Convert radians to degrees
     *
     * @method degreesToRadians
     * @param angle {Number} The angle in degrees to convert
     * @return {Number} The angle in radians
     */
    degreesToRadians: function(angle) {
        return angle * math.DEG_TO_RAD;
    },
    /**
     * Calculates the angle between two points
     *
     * @method angleBetween
     * @param pos1 {Vector|Point} The first position
     * @param pos2 {Vector|Point} The second position
     * @return {Number} The angle in radians
     */
    angleBetween: function(pos1, pos2) {
        return Math.atan2(pos2.y - pos1.y, pos2.x - pos1.x);
    }
};

module.exports = math;

},{"./random":49,"pixi.js":6}],49:[function(_dereq_,module,exports){
var support = _dereq_('../utils/support'),
    alphanum = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz';

/**
 * Static random class
 *
 * @class random
 * @static
 */
var r = {
    /**
     * This is one of the values used to generate random numbers. You should
     * never need to use this value.
     *
     * @property _x
     * @readOnly
     * @private
     * @static
     */
    _x: 0,
    /**
     * This is one of the values used to generate random numbers. You should
     * never need to use this value.
     *
     * @property _y
     * @readOnly
     * @private
     * @static
     */
    _y: 0,
    /**
     * Generates the next X value
     *
     * @method _nextX
     * @return {Number} The new X value.
     * @private
     * @static
     */
    _nextX: function() {
        return 36969 * (this._x & 0xFFFF) + (this._x >> 16);
    },
    /**
     * Generates the next Y value
     *
     * @method _nextY
     * @return {Number} The new Y value.
     * @private
     * @static
     */
    _nextY: function() {
        return 18273 * (this._y & 0xFFFF) + (this._y >> 16);
    },
    /**
     * Seeds the generator with a one, two, or three dimensional sed.
     *
     * @method seed
     * @param x {Number} The seed.
     * @param [y] {Number} The optional second dimension of the seed.
     * @param [z] {Number} The optional third dimension of the seed.
     * @static
     */
    seed: function(x, y, z) {
        this._seed = arguments;

        /* jshint -W116 */
        //use == check for undefined or null
        if(y == null && z == null) {
            this._x = x * 3253;
            this._y = this._nextX();
        } else if(z == null) {
            this._x = x * 2549 + y * 3571;
            this._y = y * 2549 + x * 3571;
        } else {
            this._x = x * 2549 + y * 3571 + z * 3253;
            this._y = x * 3253 + y * 2549 + z * 3571;
        }
        /* jshint +W116 */
    },
    /**
     * Returns the next random random value in the sequence (between 0 and 1)
     *
     * @method next
     * @alias random
     * @return {Number} The random value.
     * @static
     */
    next: function() {
        // Random number generator using George Marsaglia's MWC algorithm.

        // don't let them get stuck
        if (this._x === 0) this._x = -1;
        if (this._y === 0) this._y = -1;

        // Mix the bits.
        this._x = this._nextX();
        this._y = this._nextY();
        return ((this._x << 16) + (this._y & 0xFFFF)) / 0xFFFFFFFF + 0.5;
    },
    /**
     * Returns a random boolean based on the provided chance. The chance represents the
     * percentage chance of returning: true.
     *
     * @method bool
     * @param [chance=50] {Number} The % chance of getting true (0 - 100), defaults to 50%
     * @return {Boolean}
     * @static
     */
    bool: function(chance) {
        if(chance === undefined)
            chance = 50;

        //no chance of true
        if(chance <= 0)
            return false;

        //must always be true
        if(chance >= 100)
            return true;

        //if roll is less than change, return true
        return (r.next() * 100 < chance);
    },
    /**
     * Returns a random integer between min and max.
     *
     * @method int
     * @param [min=0] {Number} The minimun number that the result can be
     * @param [max=100] {Number} The maximun number that the result can be
     * @return {Number}
     * @static
     */
    int: function(min, max) {
        if(min !== undefined && min === max)
            return min;

        min = min || 0;
        max = max || 100;

        return Math.floor(r.next() * (max - min + 1) + min);
    },
    /**
     * Returns a random real number between min and max.
     *
     * @method real
     * @param [min=0] {Number} The minimun number that the result can be
     * @param [max=1] {Number} The maximun number that the result can be
     * @return {Number}
     * @static
     */
    real: function(min, max) {
        if(min !== undefined && min === max)
            return min;

        min = min || 0;
        max = max || 1;

        return r.next() * (max - min) + min;
    },
    /**
     * Returns a random sign based on the provided chance. The chance represents the
     * percentage chance of returning 1 (positive).
     *
     * @method sign
     * @param chance {Number} The % chance of getting positive (0 - 100), defaults to 50%
     * @return {Number} either 1 or -1
     * @static
     */
    sign: function(chance) {
        return r.bool(chance) ? 1 : -1;
    },
    /**
     * Returns a random alpha numeric string of the length specified.
     *
     * @method string
     * @param [length=16] {Number} The
     * @return {String} A random string
     * @static
     */
    string: function(length) {
        length = length || 16;
        var txt = '';

        for(var i = 0; i < length; ++i)
            txt += alphanum.charAt(Math.floor(r.next() * alphanum.length));

        return txt;
    },
    /**
     * Generates a random RFC4122 (v4) compliant UUID
     *
     * @method uuid
     * @return {String} A random uuid
     * @static
     */
    uuid: function() {
        //collect some random bytes
        var buf = r.bytes(r.__uuidBytes);

        // Per 4.4, set bits for version and `clock_seq_hi_and_reserved`
        buf[6] = (buf[6] & 0x0f) | 0x40;
        buf[8] = (buf[8] & 0x3f) | 0x80;

        var i = 0,
            bth = r.__byteToHex;

        //convert bytes to string
        return bth[buf[i++]] + bth[buf[i++]] +
                bth[buf[i++]] + bth[buf[i++]] + '-' +
                bth[buf[i++]] + bth[buf[i++]] + '-' +
                bth[buf[i++]] + bth[buf[i++]] + '-' +
                bth[buf[i++]] + bth[buf[i++]] + '-' +
                bth[buf[i++]] + bth[buf[i++]] +
                bth[buf[i++]] + bth[buf[i++]] +
                bth[buf[i++]] + bth[buf[i++]];
    },
    __uuidBytes: new Uint8Array(16),
    __byteToHex: (function() {
        var bth = [],
            htb = {};
        for (var i = 0; i < 256; i++) {
            bth[i] = (i + 0x100).toString(16).substr(1);
            htb[bth[i]] = i;
        }

        return bth;
    })(),
    /**
     * Fills a Typed Array with random bytes. If you do not pass an output param, then a default
     * Uint8Array(16) is created and returned for you.
     *
     * @method bytes
     * @param [output] {TypedArray} The output array for the random data, if none specified a new Uint8Array(16) is created
     * @static
     */
    bytes: function(ary) {
        ary = ary || (support.typedArrays ? new Uint8Array(16) : new Array(16));
        window.crypto.getRandomValues(ary);
        return ary;
    },
    /**
     * Returns a random element of an array.
     *
     * @method element
     * @param array {Array} The array to choose from
     * @param start {Number} The index of the first element to include, defaults to 0
     * @param end {Number} The index of the last element to include, defaults to array.length - 1
     * @return {Number} either 1 or -1
     * @static
     */
    element: function(array, start, end) {
        //ensure we have an array, and there are elements to check
        if(!array || !array.length)
            return null;

        //special case for 1 element
        if(array.length === 1)
            return array[0];

        //default for start
        if(!start || start < 0)
            start = start || 0;

        //default for end
        if(!end || end < 0)
            end = array.length - 1;

        return array[r.int(start, end)];
    }
};

// alias next to random
r.random = r.next;

// set the initial seed
r.seed(Math.floor(Date.now() * Math.random()));

//these polyfills are separated and exposed so that they can get tested

//if we support typed arrays we can do a good approximation of crypto.getRandomValues
r._getRandomValuesTyped = function(ary) {
    //get a Uint8 view into the buffer
    var buf = ary.buffer,
        len = buf.byteLength,
        view = new Uint8Array(buf);

    //fill the buffer one random byte at a time
    for(var i = 0, v; i < len; ++i) {
        //we only need a new random when we have pulled all the bytes out of the last one
        //which means every fourth byte we get a new random 32-bit value
        if((i & 0x03) === 0) {
            v = r.next() * 0x100000000;
        }

        //pull the next byte out of the random number
        view[i] = v >>> ((i & 0x03) << 3) & 0xff;
    }

    //return the original view which now has the data we put into the buffer
    return ary;
};

//without typed array support we can do one that returns an array of values
//but you would need to use `new Array(num)`, so there is a length
//or something like `var a = []; a[num - 1] = undefined;` so length is expanded
r._getRandomValuesArray = function(ary) {
    //fill the array with random values
    for(var i = 0; i < ary.length; ++i) {
        ary[i] = r.next() * 0x100000000;
    }

    return ary;
};

//polyfill crypto.getRandomValues if necessary
//crypto spec: http://wiki.whatwg.org/wiki/Crypto
if(!support.crypto) {
    window.crypto = window.crypto || {};

    if(support.typedArrays) {
        window.crypto.getRandomValues = r._getRandomValuesTyped;
    } else {
        window.crypto.getRandomValues = r._getRandomValuesArray;
    }
}

module.exports = r;

},{"../utils/support":70}],50:[function(_dereq_,module,exports){
var Sprite = _dereq_('../display/Sprite'),
    Texture = _dereq_('../display/Texture'),
    SpriteBatch = _dereq_('../display/SpriteBatch'),
    Vector = _dereq_('../math/Vector'),
    math = _dereq_('../math/math'),
    inherit = _dereq_('../utils/inherit'),
    C = _dereq_('../constants');

/**
 * The ParticleEmitter is the object that is placed in the world and will fire off particles based
 * on the rules and properties set on it. Generally you will want to create/use these by adding
 * them to a ParticleSystem.
 *
 * @class ParticleEmitter
 * @extends SpriteBatch
 * @constructor
 * @param name {String} The string name of the particle emitter.
 */
var ParticleEmitter = function(state, name) {
    SpriteBatch.call(this);

    this.state = state;

    /**
     * The name of the ParticleEmitter instance. This should be unique in a system, and set by the param
     * passed to the constructor.
     *
     * @property name
     * @type String
     * @readOnly
     */
    this.name = name;

    /**
     * The maximum number of particles an emitter can have active at any time.
     *
     * @property maxParticles
     * @type Number
     * @default 100
     */
    this.maxParticles = 100;

    /**
     * The width of the emitter, particles are emitted in a random integer location
     * within the width and height of the emitter.
     *
     * @property width
     * @type Number
     * @default 0
     */
    this.width = 0;

    /**
     * The height of the emitter, particles are emitted in a random integer location
     * within the width and height of the emitter.
     *
     * @property height
     * @type Number
     * @default 0
     */
    this.height = 0;

    /**
     * The default lifespan of a particle that is emitted by this ParticleEmitter, in milliseconds
     *
     * @property lifespan
     * @type Number
     * @default 2000
     */
    this.lifespan = 2000;

    /**
     * The default minSpeed of a particle that is emitted by this ParticleEmitter
     * The actual speed will be a random Vector between `minSpeed` and `maxSpeed`.
     *
     * @property minSpeed
     * @type Vector
     * @default new Vector(-100, 100)
     */
    this.minSpeed = new Vector(-100, -100);

    /**
     * The default maxSpeed of a particle that is emitted by this ParticleEmitter
     * The actual speed will be a random Vector between `minSpeed` and `maxSpeed`.
     *
     * @property maxSpeed
     * @type Vector
     * @default new Vector(100, 100)
     */
    this.maxSpeed = new Vector(100, 100);

    /**
     * The default minScale of a particle that is emitted by this ParticleEmitter
     * The actual scale will be a random number between `minScale` and `maxScale`.
     *
     * @property minScale
     * @type Number
     * @default 1
     */
    this.minScale = 1;

    /**
     * The default maxScale of a particle that is emitted by this ParticleEmitter
     * The actual scale will be a random number between `minScale` and `maxScale`.
     *
     * @property maxScale
     * @type Number
     * @default 1
     */
    this.maxScale = 1;

    /**
     * The default minRotation of a particle that is emitted by this ParticleEmitter
     * The actual rotation will be a random integer between `minRotation` and `maxRotation`.
     *
     * @property minRotation
     * @type Number
     * @default -2 * Math.PI
     */
    this.minRotation = -2 * Math.PI;

    /**
     * The default maxRotation of a particle that is emitted by this ParticleEmitter
     * The actual rotation will be a random integer between `minRotation` and `maxRotation`.
     *
     * @property maxRotation
     * @type Number
     * @default 2 * Math.PI
     */
    this.maxRotation = 2 * Math.PI;

    /**
     * The time in milliseconds between emissions of particles
     *
     * @property delay
     * @type Number
     * @default 100
     */
    this.delay = 100;

    /**
     * If true the emitter will emit particles, otherwise it will not.
     *
     * @property active
     * @type Boolean
     * @default false
     */
    this.active = false;

    this.gravity = new Vector(0, 9.87);

    //some internal trackers
    this._rate = 0; //the number of particles to emit each emission cycle
    this._total = 0; //total particles to emit

    this._emitted = 0; //total particles emitted
    this._timer = 0; //tracker for time to know when to emit particles

    //params for particle ctor
    this._particle = null; //the Sprite object to use as the "template" particle
    this._textures = null; //the textures to choose from when getting a particle
    this._pool = []; //the pool to release particles into when they are done
};

inherit(ParticleEmitter, SpriteBatch, {
    /**
     * Starts the particle emission, must call `setup` first to setup
     * what kind of particle to emit.
     *
     * @method start
     * @param [lifespan=2000] {Number} The lifespan of a particle in ms
     * @param [rate=1] {Number} The number of particles to emit each emission, Infinity means to operate in "burst" mode.
     *      Burst mode is when the emitter dumps all the particles at once, then goes inactive. The "delay" param is
     *      ignored when the emitter is in burst mode.
     * @param [total=100] {Number} The total number of particles to emit
     * @param [delay=250] {Number} The time between each particle emission in ms
     * @return {ParticleEmitter} Returns itself.
     * @chainable
     */
    start: function(lifespan, rate, total, delay) {
        this.active = true;

        this.lifespan = lifespan || 2000;
        this.delay = delay || 250;
        this._rate = rate || 1;
        this._total = total || 100;

        this._timer = 0;

        //special case for burst mode
        if(this._rate === Infinity) {
            this.active = false;

            for(var i = 0; i < this._total; ++i) {
                this.emitParticle();
            }
        }

        return this;
    },
    /**
     * Deactivates the emitter. Particles that are already emitted will continue to
     * decay and die, but no new particles will be emitted.
     *
     * @method stop
     * @return {ParticleEmitter} Returns itself.
     * @chainable
     */
    stop: function() {
        this.active = false;

        return this;
    },
    /**
     * Sets up the particles to be emitted
     *
     * @method setup
     * @param sprite {Sprite|Array<Texture>|Texture} Pass a sprite to be clones as a particle,
     *      or an array of textures to be randomly chosen from for different particles,
     *      or a single texture to use for each particle.
     * @param [collide=gf.DIRECTION.ALL] {Number} The directions the particles are allowed to collide in, use gf.DIRECTION bit flags
     * @return {ParticleEmitter} Returns itself.
     * @chainable
     */
    setup: function(sprite, collide) {
        if(collide === undefined)
            collide = C.DIRECTION.ALL;

        //single texture
        if(sprite instanceof Texture) {
            this._particle = new Sprite(sprite);
            this._textures = [sprite];
        }
        //array of textures
        else if(Array.isArray(sprite)) {
            this._particle = new Sprite(sprite[0]);
            this._textures = sprite;
        }
        //an actual sprite
        else {
            this._particle = sprite;
            this._textures = [sprite.texture];
        }

        return this;
    },
    /**
     * Gets a particle from the pool and sets it up.
     *
     * @method _get
     * @return {Sprite} The particle to use
     * @private
     */
    _get: function() {
        if(this._emitted >= this._total || this._emitted > this.maxParticles)
            return null;

        var spr = this._pool.pop();

        if(!spr) {
            spr = this._particle.clone();
        }

        spr.setTexture(math.rand.element(this._textures));
        spr.visible = true;

        this.addChild(spr);
        this._emitted++;

        return spr;
    },
    /**
     * Frees a particle back into the pool and hides it.
     *
     * @method _free
     * @param sprite {Sprite} The particle to free
     * @private
     */
    _free: function(spr) {
        spr.visible = false;
        this._pool.push(spr);

        this.removeChild(spr);
        this._emitted--;
    },
    /**
     * Emits a single particle and sets the position, scale, lifespan, and velocity
     *
     * @method emitParticle
     * @return {ParticleEmitter} Returns itself.
     * @chainable
     */
    emitParticle: function() {
        var part = this._get();

        if(!part)
            return;

        //set optionally random position
        part.setPosition(
            math.rand.int(0, this.width),
            math.rand.int(0, this.height)
        );

        //set scale
        part.scale.x = part.scale.y = math.rand.real(this.minScale, this.maxScale);

        //set lifespan
        part.lifespan = this.lifespan;

        //set velocity
        part.setVelocity(
            math.rand.int(this.minSpeed.x, this.maxSpeed.x),
            math.rand.int(this.minSpeed.y, this.maxSpeed.y)
        );

        //part.body.angularVelocity = math.rand.int(this.minRotation, this.maxRotation);

        return this;
    },
    /**
     * Called internally by the ParticleSystem each frame to update each particle's lifespan.
     *
     * @method update
     * @param dt {Number} The number of seconds that have passed since last call
     * @private
     */
    update: function(dt) {
        var t = dt * 1000;

        //update each of the particle lifetimes
        for(var c = 0; c < this.children.length; ++c) {
            var child = this.children[c];
            child.lifespan -= t;

            child.position.x += (child._velocity.x += this.gravity.x) * dt;
            child.position.y += (child._velocity.y += this.gravity.y) * dt;

            if(child.lifespan <= 0) {
                this.emit('particle.expire', child);
                this._free(child);
            }
        }

        //if no longer active, we are done here
        if(!this.active)
            return;

        //increment time we have waited
        this._timer += t;

        //if we waited more than delay, emit some particles
        if(this._timer >= this.delay) {
            this._timer -= this.delay;

            for(var i = 0; i < this._rate; ++i) {
                this.emitParticle();
            }
        }
    }
});

module.exports = ParticleEmitter;

},{"../constants":11,"../display/Sprite":19,"../display/SpriteBatch":20,"../display/Texture":21,"../math/Vector":47,"../math/math":48,"../utils/inherit":69}],51:[function(_dereq_,module,exports){
var Emitter = _dereq_('./ParticleEmitter'),
    Container = _dereq_('../display/Container'),
    inherit = _dereq_('../utils/inherit');

/**
 * The ParticleSystem controls the system of particle emitters and their particles. It contains all the emitters
 * and updates them each frame. An instance of this is created for you in a world instance, which is a property
 * of a game state. The general usage for this class is:
 *
 * @class ParticleSystem
 * @extends Container
 * @constructor
 */
var ParticleSystem = function(state) {
    Container.call(this);

    this.state = state;

    /**
     * The emitters that are contained in this system, keyed by name
     *
     * @property emitter
     * @type Object
     * @readOnly
     */
    this.emitters = {};

    /**
     * The next ID to use for an emitter where no name was passed
     *
     * @property nextId
     * @type Number
     * @private
     */
    this.nextId = 0;
};

inherit(ParticleSystem, Container, {
    /**
     * Adds a particle emitter to the system, creating one if necessary.
     *
     * There are 3 ways to use this function to add an emitter to the system. The simplest case
     * is to pass a string for the name, and let the manager create a normal gf.ParticleEmitter for you
     * with the name you provided. The second usage is to pass a class that is a decendant of gf.ParticleEmitter.
     *
     * For example:
     *
     * ```
     * function MyEmitter(name) {
     *     gf.ParticleEmitter.call(name);
     * }
     * gf.inherit(MyEmitter, gf.ParticleEmitter);
     *
     * game.world.particles.add(MyEmitter); //adds a new instance of your emitter
     * ```
     *
     * The final usage is to pass an Emitter that is already created. In this case the system will
     * add the emitter to the list based on `emitter.name`.
     *
     * @method add
     * @param emitter {String|Function|ParticleEmitter} The emitter name, constructor, or emitter instance to add.
     * @return {ParticleEmitter} The emitter that was added
     */
    add: function(Name) {
        var emitter;

        //create an emitter if a string is passed
        if(typeof Name === 'string') {
            emitter = new Emitter(this.state, Name);
        }
        //create an emitter of the instance passed
        else if(typeof Name === 'function') {
            emitter = new Name(this.state);
        }
        //a pre-created emitter, ensure game is set correctly
        else {
            emitter = Name;
        }

        if(!emitter.name)
            emitter.name = 'emitter_' + (this.nextId++);

        this.emitters[emitter.name] = emitter;
        this.addChild(emitter);

        return emitter;
    },
    /**
     * Removes an emitter from the system
     *
     * @method remove
     * @param emitter {String|ParticleEmitter} The name of the emitter to remove, or the emitter instance itself.
     * @return {ParticleSystem} Returns itself.
     * @chainable
     */
    remove: function(emitter) {
        if(typeof emitter === 'string')
            emitter = this.emitters[emitter];

        if(emitter.parent)
            emitter.parent.removeChild(emitter);

        delete this.emitters[emitter.name];

        return this;
    },
    /**
     * Called internally by the World each frame to update each Particle Emitter
     *
     * @method update
     * @param dt {Number} The number of seconds that have passed since last call
     * @private
     */
    updateTransform: function() {
        Container.prototype.updateTransform.apply(this, arguments);

        //get delta
        var dt = this.state.game.timings.lastDelta;

        for(var i = 0, il = this.children.length; i < il; ++i) {
            var emitter = this.children[i];

            if(emitter.update)
                emitter.update(dt);
        }
    }
});

module.exports = ParticleSystem;

},{"../display/Container":16,"../utils/inherit":69,"./ParticleEmitter":50}],52:[function(_dereq_,module,exports){
var Rectangle = _dereq_('../geom/Rectangle'),
    Circle = _dereq_('../geom/Circle'),
    Polygon = _dereq_('../geom/Polygon'),
    Vector = _dereq_('../math/Vector'),
    Tile = _dereq_('../tilemap/Tile'),
    math = _dereq_('../math/math'),
    inherit = _dereq_('../utils/inherit'),
    cp = _dereq_('chipmunk');

/**
 * The PhysicsSystem is the wrapper around the chipmunk-js physics library that integrates
 * grapefruit objects into the physics world. It is in charge of managing objects in the physics
 * space. Generally you would not create this yourself and instead would use the `.physics` property
 * of a State.
 *
 * @class PhysicsSystem
 * @extends Object
 * @constructor
 * @param state {State} The state instance this system belongs to.
 * @param [options] {Object} The options for the physics system.
 * @param [options.gravity=new Vector(0, 9.87)] {Vector} The gravity of the space
 */
var PhysicsSystem = function(state, options) {
    //default options
    options = options || {};
    options.gravity = options.gravity instanceof Vector ? options.gravity : new Vector(0, 9.87);
    options.iterations = options.iterations || 10;
    options.sleepTimeThreshold = options.sleepTimeThreshold !== undefined ? options.sleepTimeThreshold : 0.5;
    options.collisionSlop = options.collisionSlop !== undefined ? options.collisionSlop : 0.1;
    options.stepTime = options.stepTime || (1 / 60);

    /**
     * The state instance this system belongs to
     *
     * @property state
     * @type State
     */
    this.state = state;

    /**
     * The delta time to use as the constant for physics simulation
     *
     * @property stepTime
     * @type Number
     * @default (1 / 60)
     */
    this.stepTime = options.stepTime;

    /**
     * The chipmunk space instance that will run all the physics simulations
     *
     * @property space
     * @type cp.Space
     * @readOnly
     */
    this.space = new cp.Space();

    /**
     * The gravity of the physics space
     *
     * @property gravity
     * @type Vector
     */
    this.gravity = this.space.gravity = options.gravity;

    //Time a body must remain idle to fall asleep
    //see: http://chipmunk-physics.net/release/ChipmunkLatest-API-Reference/structcp_space.html#a928d74741904aae266a9efff5b5f68f7
    this.space.sleepTimeThreshold = options.sleepTimeThreshold;

    //Amount of encouraged penetration between colliding shapes.
    //see: http://chipmunk-physics.net/release/ChipmunkLatest-API-Reference/structcp_space.html#af1bec644a24e12bfc642a942a88520f7
    this.space.collisionSlop = options.collisionSlop;

    //These two collision scenarios are separate because we don't
    //want tiles to collide with tiles all the time

    //sprite - sprite collisions
    this.space.addCollisionHandler(
        PhysicsSystem.COLLISION_TYPE.SPRITE,
        PhysicsSystem.COLLISION_TYPE.SPRITE,
        this.onCollisionBegin.bind(this), //begin
        null, //preSolve
        this.onCollisionPostSolve.bind(this), //postSolve
        this.onCollisionEnd.bind(this) //separate
    );

    //sprite - tile collisions
    this.space.addCollisionHandler(
        PhysicsSystem.COLLISION_TYPE.SPRITE,
        PhysicsSystem.COLLISION_TYPE.TILE,
        this.onCollisionBegin.bind(this), //begin
        null, //preSolve
        this.onCollisionPostSolve.bind(this), //postSolve
        this.onCollisionEnd.bind(this) //separate
    );

    /**
     * The actions to perform at the next post-step callback. These are
     * for actions (like add, remove) that cannot be performed during a
     * simulation step.
     *
     * @property actonQueue
     * @type Array<Object>
     * @private
     */
    this.actionQueue = [];

    /**
     * The callback functions to call on the next frame
     *
     * @property tickCallbacks
     * @type Array<Function>
     * @private
     */
    this.tickCallbacks = [];

    /**
     * The number of steps to skip, tracks `this.skip(num)`
     *
     * @property _skip
     * @type Number
     * @private
     */
    this._skip = 0;

    this._updateNum = 0;
    this._paused = false;
};

inherit(PhysicsSystem, Object, {
    /**
     * Pauses physics simulation
     *
     * @method pause
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    pause: function() {
        this._paused = true;

        return this;
    },
    /**
     * Resumes physics simulation
     *
     * @method resume
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    resume: function() {
        this._paused = false;

        return this;
    },
    /**
     * Skips the specified number of frame steps
     *
     * @method skip
     * @param num {Number} Number of steps to skip
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    skip: function(num) {
        this._skip += num;

        return this;
    },
    /**
     * Skips the next frame step
     *
     * @method skipNext
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    skipNext: function() {
        return this.skip(1);
    },
    /**
     * Registers a callback to be executed on the next frame step
     *
     * @method nextTick
     * @param fn {Function} The callback to register
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     * @async
     */
    nextTick: function(fn) {
        this.tickCallbacks.push(fn);

        return this;
    },
    /**
     * Returns the collision type of a sprite
     *
     * @method getCollisionType
     * @param spr {Sprite} The sprite to check
     * @return {Number} The collision type
     */
    getCollisionType: function(spr) {
        if(spr instanceof Tile) {
            return PhysicsSystem.COLLISION_TYPE.TILE;
        } else {
            return PhysicsSystem.COLLISION_TYPE.SPRITE;
        }
    },
    /**
     * Adds a sprite to the physics simulation
     *
     * @method add
     * @param spr {Sprite} The sprite to add
     * @param [callback] {Function} The callback to call once the sprite has been added
     * @return {Sprite} The sprite that was added
     * @async
     */
    add: function(spr, cb) {
        //already in space with body(s)
        if(spr._phys.active)
            return;

        var body = this._createBody(spr),
            shape = this._createShape(spr, body);

        spr._phys.active = true;
        this.actionQueue.push(['add', {
            spr: spr,
            body: body,
            shape: shape
        }, cb]);
        this.act();

        return spr;
    },
    addControlBody: function(spr, cb) {
        //see Chipmunk2D Tank Demo:
        //https://github.com/slembcke/Chipmunk2D/blob/master/Demo/Tank.c#L106
        var body = spr._phys.body;

        if(!body.isStatic()) {
            var cbody = new cp.Body(Infinity, Infinity),
                cpivot = new cp.PivotJoint(cbody, body, cp.vzero, cp.vzero),
                cgear;

            cpivot.maxBias = 0; //disable join correction
            cpivot.maxForce = 10000; //emulate linear friction

            //infinite inertia cannot rotate, so we wouldn't need a gear joint
            if(body.i !== Infinity) {
                cgear = new cp.GearJoint(cbody, body, 0, 1);
                cgear.errorBias = 0; //attempt to fully correct the joint each step
                cgear.maxBias = 1.2; //but limit the angular correction
                cgear.maxForce = 50000; //emulate angular friction
            }

            this.actionQueue.push(['addControl', {
                spr: spr,
                body: cbody,
                pivot: cpivot,
                gear: cgear
            }, cb]);
            this.act();
        }

        return spr;
    },
    /**
     * Removes a sprite from the physics simulation
     *
     * @method remove
     * @param spr {Sprite} The sprite to remove
     * @param [callback] {Function} The callback to call once the sprite has been removed
     * @return {Sprite} The sprite that was removed
     * @async
     */
    remove: function(spr, cb) {
        if(!spr || !spr._phys.active)
            return;

        spr._phys.active = false;
        this.actionQueue.push(['remove', spr._phys, cb]);
        this.act();

        return spr;
    },
    /**
     * Reindexes a sprite's shape in the simulation, useful if it looks
     * like changes are being cached.
     *
     * @method reindex
     * @param spr {Sprite} The sprite to reindex
     * @param [callback] {Function} The callback to call once the sprite has been reindexed
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     * @async
     */
    reindex: function(spr, cb) {
        if(!spr || !spr._phys.active)
            return;

        this.actionQueue.push(['reindex', spr._phys.shape, cb]);
        this.act();

        return this;
    },
    /**
     * Reindexes all static bodies in the simulation.
     *
     * @method reindexStatic
     * @param [callback] {Function} The callback to call once reindexing completes
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     * @async
     */
    reindexStatic: function(cb) {
        this.actionQueue.push(['reindexStatic', null, cb]);
        this.act();

        return this;
    },
    /**
     * Adds a custom shape to a sprite, useful for a single sprite to have multiple
     * different collision shapes (including sensors).
     *
     * @method addCustomShape
     * @param spr {Sprite} The sprite to add the shape to
     * @param poly {Circle|Rectangle|Polygon} The shape to create
     * @param sensor {Boolean} Is this a sensor shape, if so you will get a collision callback, but no solve
     * @param [callback] {Function} The callback to call once the shape has been added
     * @return {cp.Shape} The shape that was created
     * @async
     */
    addCustomShape: function(spr, poly, sensor, cb) {
        if(!spr || !spr._phys.body)
            return;

        var s = this._createShape(spr, spr._phys.body, poly);

        s.width = spr.width;
        s.height = spr.height;
        s.sprite = spr;

        s.setSensor(sensor);
        s.setElasticity(0);
        s.setSensor(sensor !== undefined ? sensor : spr.sensor);
        s.setCollisionType(this.getCollisionType(spr));
        s.setFriction(spr.friction || 0);

        this.actionQueue.push(['addCustomShape', { spr: spr, shape: s }, cb]);
        this.act();

        return s;
    },
    /**
     * Sets the mass of a sprite's physics body.
     *
     * @method setMass
     * @param spr {Sprite} The sprite to set the mass for
     * @param mass {Number} The mass to set
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    setMass: function(spr, mass) {
        if(!spr || !spr._phys.body)
            return;

        spr._phys.body.setMass(mass);

        return this;
    },
    /**
     * Sets the velocity of a sprite's physics body.
     *
     * @method setVelocity
     * @param spr {Sprite} The sprite to set the velocity for
     * @param velocity {Vector} The velocity to set to
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    setVelocity: function(spr, vel) {
        if(!spr)
            return;

        //update control body velocity (and pivot contraint makes regular follow)
        if(spr._phys.control) {
            spr._phys.control.body.setVel(vel);
        }
        //if no control body then update real body
        else if(spr._phys.body) {
            spr._phys.body.setVel(vel);
        }

        return this;
    },
    /**
     * Sets the position of a sprite's physics body.
     *
     * @method setPosition
     * @param spr {Sprite} The sprite to set the position for
     * @param position {Vector} The position to set to
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    setPosition: function(spr, pos) {
        if(!spr)
            return;

        //update body position
        if(spr._phys.body) {
            spr._phys.body.setPos(pos.clone());
        }

        //update control body position
        if(spr._phys.control) {
            spr._phys.control.body.setPos(pos.clone());
        }


        return this;
    },
    /**
     * Sets the rotation of a sprite's physics body.
     *
     * @method setRotation
     * @param spr {Sprite} The sprite to set the rotation for
     * @param rotation {Number} The rotation to set to in radians
     * @return {PhysicsSystem} Returns itself.
     * @chainable
     */
    setRotation: function(spr, rads) {
        if(!spr)
            return;

        //update control body rotation (and gear contraint makes regular follow)
        if(spr._phys.control) {
            spr._phys.control.body.setAngle(rads);
        }
        //if no control body then update real body
        else if(spr._phys.body) {
            spr._phys.body.setAngle(rads);
        }

        return this;
    },
    /**
     * Called each physics step iteration. This is detached from the frame rendering and
     * runs at a constant step.
     *
     * @method update
     * @private
     */
    update: function(dt) {
        if(this._paused)
            return;

        while(this.tickCallbacks.length)
            (this.tickCallbacks.shift()).call(this);

        if(this._skip)
            return this._skip--;

        //execute the physics step
        this.space.step(this.stepTime);

        //go through each changed shape
        var alpha = dt / this.stepTime,
            num = this._updateNum++,
            spr, body;

        if(alpha > 1) {
            this.update(dt - this.stepTime);
            alpha = 1;
        }

        this.space.activeShapes.each(function(shape) {
            body = shape.body;
            spr = shape.sprite;

            //already updated this body
            if(body._updateNum === num)
                return;

            body._updateNum = num;

            //update sprite
            spr.position.lerp(body.p, alpha).round();
            spr.rotation += (body.a - spr.rotation) * alpha;
            spr.rotation = math.round(spr.rotation);

            //the sprite has changed due to a physics update, emit that event
            spr.emit('physUpdate');
        });
    },
    /**
     * Called when a collision begins in the system
     *
     * @method onCollisionBegin
     * @param arbiter {cp.Arbiter} The arbiter of the collision
     * @param space {cp.Space} The space the collision occurred in
     * @private
     */
    onCollisionBegin: function(arbiter) {//, space) {
        var shapes = arbiter.getShapes(),
            spr1 = shapes[0].sprite,
            spr2 = shapes[1].sprite;

        //only call the sensor collisions here
        if(shapes[0].sensor || shapes[1].sensor) {
            spr1.onCollision(spr2, arbiter.getNormal(0), shapes[1], shapes[0]);
            spr2.onCollision(spr1, arbiter.getNormal(0), shapes[0], shapes[1]);
        }

        //maintain the colliding state
        return true;
    },
    /**
     * Called after a collision is solved in the system
     *
     * @method onCollisionPostSolve
     * @param arbiter {cp.Arbiter} The arbiter of the collision
     * @param space {cp.Space} The space the collision occurred in
     * @private
     */
    onCollisionPostSolve: function(arbiter) {//, space) {
        var shapes = arbiter.getShapes(),
            spr1 = shapes[0].sprite,
            spr2 = shapes[1].sprite;

        if(arbiter.isFirstContact()) {
            spr1.onCollision(spr2, arbiter.totalImpulse(), shapes[1], shapes[0]);
            spr2.onCollision(spr1, arbiter.totalImpulse(), shapes[0], shapes[1]);
        }

        //maintain the colliding state
        return true;
    },
    /**
     * Called after a collision ends in the system (separation)
     *
     * @method onCollisionEnd
     * @param arbiter {cp.Arbiter} The arbiter of the collision
     * @param space {cp.Space} The space the collision occurred in
     * @private
     */
    onCollisionEnd: function(arbiter) {//, space) {
        var shapes = arbiter.getShapes(),
            spr1 = shapes[0].sprite,
            spr2 = shapes[1].sprite;

        spr1.onSeparate(spr2, shapes[1], shapes[0]);
        spr2.onSeparate(spr1, shapes[0], shapes[1]);

        //maintain the colliding state
        return true;
    },
    /**
     * Attempts to perform the postStep actions that have been queued. If the space
     * is currently locked, then it waits until after the step to run the actions.
     *
     * @method onCollisionEnd
     * @private
     */
    act: function() {
        if(this.space.locked) {
            this.space.addPostStepCallback(this.onPostStep.bind(this));
        } else {
            //for async behavior
            var self = this;
            setTimeout(function() {
                self.onPostStep();
            }, 1);
        }
    },
    /**
     * Processes the action queue after a step is unlocked.
     *
     * @method onPostStep
     * @private
     */
    onPostStep: function() {
        //remove items
        while(this.actionQueue.length) {
            var a = this.actionQueue.shift(),
                act = a[0],
                data = a[1],
                cb = a[2];

            switch(act) {
                case 'add':
                    data.body.setPos(data.spr.position.clone());
                    if(!data.body.isStatic()) {
                        this.space.addBody(data.body);
                    }

                    this.space.addShape(data.shape);

                    data.spr._phys.body = data.body;
                    data.spr._phys.shape = data.shape;
                    data.body.sprite = data.spr;
                    break;

                case 'addControl':
                    data.body.setPos(data.spr.position.clone());
                    this.space.addBody(data.body);
                    this.space.addConstraint(data.pivot);

                    if(data.gear) this.space.addConstraint(data.gear);

                    data.spr._phys.control = data;
                    delete data.spr; //no need for that extra reference to lay around
                    break;

                case 'addCustomShape':
                    if(!data.spr._phys.customShapes) {
                        data.spr._phys.customShapes = [];
                    }

                    data.spr._phys.customShapes.push(data.shape);
                    this.space.addShape(data.shape);
                    break;

                case 'remove':
                    if(data.body.space) {
                        this.space.removeBody(data.body);
                    }

                    if(data.shape.space) {
                        this.space.removeShape(data.shape);
                    }

                    if(data.control) {
                        if(data.control.body.space) {
                            this.space.removeBody(data.control.body);
                        }

                        if(data.control.pivot.space) {
                            this.space.removeConstraint(data.control.pivot);
                        }

                        if(data.control.gear && data.control.gear.space) {
                            this.space.removeConstraint(data.control.gear);
                        }
                    }

                    if(data.customShapes) {
                        for(var i = data.customShapes.length - 1; i > -1; --i) {
                            this.space.removeShape(data.customShapes[i]);
                        }
                    }

                    //remove references
                    data.body = null;
                    data.shape.sprite = null;
                    data.shape = null;
                    data.customShapes = null;
                    break;

                case 'reindex':
                    this.space.reindexShape(data);
                    break;

                case 'reindexStatic':
                    this.space.reindexStatic();
                    break;
            }

            if(cb)
                cb.call(this);
        }
    },
    /**
     * Creates a physics body for a sprite
     *
     * @method _createBody
     * @param spr {Sprite} The sprite to create a body for
     * @return {cp.Body} The chipmunk-js physics body
     * @private
     */
    _createBody: function(spr) {
        var mass = spr.mass || 1,
            inertia = spr.inertia || cp.momentForBox(mass, spr.width, spr.height) || Infinity,
            body = new cp.Body(mass, inertia);

        if(mass === Infinity) {
            body.nodeIdleTime = Infinity;
        } else {
            body.nodeIdleTime = 0;
        }

        return body;
    },
    /**
     * Creates a collision shape for a sprite
     *
     * @method _createShape
     * @param spr {Sprite} The sprite to create a shape for
     * @param body {cp.Body} The body to attach the shape to
     * @param [poly] {Circle|Rectangle|Polygon} The shape to create, defaults to `spr.hitArea`
     * @return {cp.Shape} The chipmunk-js collision shape
     * @private
     */
    _createShape: function(spr, body, poly) {
        var shape,
            hit = poly || spr.hitArea,
            ax = spr.anchor ? spr.anchor.x : 0,
            ay = spr.anchor ? spr.anchor.y : 0,
            aw = spr.width * ax,
            ah = spr.height * ay;

        //specified shape
        if(hit) {
            if(hit instanceof Rectangle) {
                //convert the top-left anchored rectangle to left,right,bottom,top values
                //for chipmunk space that will corospond to our own
                var l = hit.x,
                    r = hit.x + hit.width,
                    b = hit.y - spr.height,
                    t = b + hit.height;

                l -= aw;
                r -= aw;

                b += spr.height - ah;
                t += spr.height - ah;

                shape = new cp.BoxShape2(body, new cp.BB(l, b, r, t));
            }
            else if(hit instanceof Circle) {
                //the offset needs to move the circle to the sprite center based on the sprite's anchor (bottom-left)
                var offset = new Vector(
                    ((spr.width / 2) - aw) + hit.x,
                    ((spr.height / 2) - ah) + hit.y
                );

                shape = new cp.CircleShape(body, hit.radius, offset);
            }
            else if(hit instanceof Polygon) {
                //cp shapes anchors are 0.5,0.5, but a polygon uses 0,0 as the topleft
                //of the bounding rect so we have to convert
                var points = [],
                    ps = hit.points;

                for(var i = 0; i < ps.length; ++i) {
                    var p = ps[i];

                    points.push(p.x - aw);
                    points.push(p.y - ah);
                }

                shape = new cp.PolyShape(body, cp.convexHull(points, null, 0), cp.vzero);
            }
        }

        //default box shape
        if(!shape) {
            shape = new cp.BoxShape2(body, new cp.BB(0, -spr.height, spr.width, 0));
        }

        //this.space.addShape(shape);

        shape.width = spr.width;
        shape.height = spr.height;
        shape.sprite = spr;
        shape.setElasticity(spr.bounce || spr.elasticity || 0);
        shape.setSensor(spr.sensor);
        shape.setCollisionType(this.getCollisionType(spr));
        shape.setFriction(spr.friction || 0);

        shape.group = spr.shapeGroup || 0;

        return shape;
    }
});

PhysicsSystem.COLLISION_TYPE = {
    SPRITE: 0,
    TILE: 1
};

module.exports = PhysicsSystem;

},{"../geom/Circle":34,"../geom/Polygon":36,"../geom/Rectangle":37,"../math/Vector":47,"../math/math":48,"../tilemap/Tile":57,"../utils/inherit":69,"chipmunk":1}],53:[function(_dereq_,module,exports){
var Vector = _dereq_('../math/Vector');

/**
 * Physics mixin. This will add physics capabilities to the class it mixes into.
 *
 * @class PhysicsTarget
 * @constructor
 */
module.exports = function() {
    /**
     * The physics namespace that all physics properties go into. Those properties are:
     *  - system {PhysicsSystem} PhysicsSystem that this object is a part of.
     *  - active {Boolean} Whether or not this target is actively having physics simulated.
     *
     * @property _phys
     * @type Object
     * @default {}
     * @private
     * @readOnly
     */
    this._phys = {};

    /**
     * The mass of this object, please use setMass to set this value
     *
     * @property mass
     * @type Number
     * @default 0
     * @readOnly
     */
    this.mass = 0;

    /**
     * The moment of inertia of this object, only set this before enabling physics (has no effect after enabling)
     *
     * @property inertia
     * @type Number
     * @default 0
     */
    this.inertia = 0;

    /**
     * The internal velocity value, used as a reusable vector for the setVelocity function. Setting
     * this directly *has no effect*.
     *
     * @property _velocity
     * @type Vector
     * @private
     * @readOnly
     */
    this._velocity = new Vector();

    /**
     * Enables physics for this sprite
     *
     * @method enablePhysics
     * @param system {PhysicsSystem} The system for the sprite to be in
     * @param callback {Function} The callback function to call after physics have been enabled
     * @return {mixed} Returns itself.
     * @chainable
     * @async
     */
    this.enablePhysics = function(sys, cb) {
        var self = this;

        if(typeof sys === 'function') {
            cb = sys;
            sys = null;
        }

        //is a system is passed use it
        if(sys) {
            //if active, remove from current system
            if(this._phys.active) {
                //remove from old system
                this._phys.system.remove(this, function() {
                    //add to new system when completed
                    sys.add(self, cb);
                });
            }
            //if inactive add to new system immediately
            else {
                sys.add(this, cb);
            }

            //reassign new system
            this._phys.system = sys;
        }
        //if no system passed (or same one passed) just add to current stored system
        else {
            this._phys.system.add(this, cb);
        }

        return this;
    };

    /**
     * Disbales physics for this sprite
     *
     * @method disablePhysics
     * @param callback {Function} The callback function to call after physics have been disabled
     * @return {mixed} Returns itself.
     * @chainable
     * @async
     */
    this.disablePhysics = function(cb) {
        //if we have a cached system, remove from it
        if(this._phys.system) {
            this._phys.system.remove(this, cb);
        }

        return this;
    };

    /**
     * Reindexes the collisions for this sprite, useful when moving the sprite a great deal
     * (like to a new world)
     *
     * @method reindex
     * @param callback {Function} The callback function to call after the sprite has been reindexed
     * @return {mixed} Returns itself.
     * @chainable
     * @async
     */
    this.reindex = function(cb) {
        //if we have a cached system, reindex
        if(this._phys.system) {
            this._phys.system.reindex(this, cb);
        }

        return this;
    };

    /**
     * Sets the mass of this sprite
     *
     * @method setMass
     * @param mass {Number} The new mass of the object
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.setMass = function(mass) {
        if(this._phys.system) {
            this._phys.system.setMass(this, mass);
        }

        return this;
    };

    /**
     * Sets the velocity of this sprite
     *
     * @method setVelocity
     * @param x {Number|Vector} The x coord to set the velocity to, if a Vector is passed the y param is ignored
     * @param y {Number} The y coord to set the velocity to
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.setVelocity = function(x, y) {
        y = x.y !== undefined ? x.y : (y || 0);
        x = x.x !== undefined ? x.x : (x || 0);

        this._velocity.set(x, y);

        if(this._phys.system) {
            this._phys.system.setVelocity(this, this._velocity);
        }

        return this;
    };

    /**
     * Sets the rotation of this sprite
     *
     * @method setRotation
     * @param rotation {Number} The new rotation of the object in radians
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.setRotation = function(rads) {
        this.rotation = rads;

        if(this._phys.system) {
            this._phys.system.setRotation(this, rads);
        }

        return this;
    };

    /**
     * Sets the position of this sprite
     *
     * @method setPosition
     * @param x {Number|Vector} The x coord to set position to, if a Vector is passed the y param is ignored
     * @param y {Number} The y coord to set position to
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.setPosition = function(x, y) {
        y = x.y !== undefined ? x.y : (y || 0);
        x = x.x !== undefined ? x.x : (x || 0);

        this.position.set(x, y);

        if(this._phys.system) {
            this._phys.system.setPosition(this, this.position);
        }

        return this;
    };

    /**
     * On Collision Event
     *      called when this sprite collides into another, or is being collided into by another.
     *      By default if something collides with a collectable sprite we destroy the collectable
     *      and if we collide with a solid tile we kill our velocity. This method will emit a
     *      'collision' event that you can listen for
     *
     * @event collision
     * @param obj {Sprite} Colliding sprite
     * @param vec {Vector} Collision vector (for sensors this is normalized)
     * @param colShape {cp.Shape} The colliding physics shape
     * @param myShape {cp.Shape} Your physics shape that caused the collision
     */
    this.onCollision = function(obj, vec, colShape, myShape) {
        this.emit('collision', obj, vec, colShape, myShape);
    };

    /**
     * On Seperate Event
     *      called when this sprite collides into another, or is being collided into by another.
     *      By default if something collides with a collectable sprite we destroy the collectable
     *      and if we collide with a solid tile we kill our velocity. This method will emit a
     *      'collision' event that you can listen for
     *
     * @event separate
     * @param obj {Sprite} Colliding sprite
     * @param colShape {cp.Shape} The colliding physics shape
     * @param myShape {cp.Shape} Your physics shape that caused the collision
     */
    this.onSeparate = function(obj, colShape, myShape) {
        this.emit('separate', obj, colShape, myShape);
    };
};

},{"../math/Vector":47}],54:[function(_dereq_,module,exports){
var SpriteBatch = _dereq_('../display/SpriteBatch'),
    ObjectPool = _dereq_('../utils/ObjectPool'),
    Texture = _dereq_('../display/Texture'),
    Sprite = _dereq_('../display/Sprite'),
    Vector = _dereq_('../math/Vector'),
    Rectangle = _dereq_('../geom/Rectangle'),
    utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    PIXI = _dereq_('pixi.js');

/**
 * A Text Object will create (a) line(s) of text using bitmap font. To split a line you can use "\n", "\r" or "\r\n"
 * You can generate the fnt files using [bmfont](http://www.angelcode.com/products/bmfont/) for windows or
 * [bmglyph](http://www.bmglyph.com/) for mac.
 *
 * @class BitmapText
 * @extends SpriteBatch
 * @constructor
 * @param text {String} The copy that you would like the text to display
 * @param font {Object} The font data object (this is generally grabbed from `game.cache.getBitmapFont('mykey')`);
 * @param font.name {String} The name of the font
 * @param font.size {Number} The base size of the font
 * @param font.lineHeight {Number} The line height of the font
 * @param font.chars {Object} The characters in the font, each should have a texture and kerning info
 * @param [style] {Object} The style parameters
 * @param [style.size=null] {Number} The font size of the text, overrides the font's size
 * @param [style.align="left"] {String} An alignment of the multiline text ("left", "center" or "right")
 */
var BitmapText = function(text, font, style) {
    SpriteBatch.call(this);

    /**
     * Whether or not the bitmap text is dirty and should be redrawn
     *
     * @property dirty
     * @type Boolean
     * @default false
     */
    this.dirty = false;

    /**
     * The font descriptor that holds the font data (name, size, chars, etc)
     *
     * @property font
     * @type Object
     * @readOnly
     */
    this.font = font;

    /**
     * A monospacing to apply to the font instead of the actual character/kerning spacing.
     * When set to `0` the default font values will be used.
     *
     * @property monospace
     * @type Number
     * @default 0
     */
    this.monospace = 0;

    /**
     * The actual text that is currently rendered, please use the `text` property
     * and do not set this directly.
     *
     * @property _text
     * @type String
     * @readOnly
     * @private
     */
    this._text = text;

    /**
     * The current style of the bitmap text
     *
     * @property _style
     * @type Object
     * @readOnly
     * @private
     */
    this._style = {
        size: null,
        align: 'left'
    };

    /**
     * The sprite pool to grab character sprites from
     *
     * @property sprites
     * @type ObjectPool
     * @readOnly
     * @private
     */
    this.sprites = new ObjectPool(Sprite, this);

    this.text = text;
    this.style = style;
};

inherit(BitmapText, SpriteBatch, {
    /**
     * Renders the text character sprites when the text is dirty. This is
     * automatically called when the text/style becomes dirty.
     *
     * @method renderText
     */
    renderText: function() {
        var font = this.font,
            pos = new Vector(),
            prevCode = null,
            chars = [],
            maxLineWidth = 0,
            lineWidths = [],
            line = 0,
            scale = this._style.size / font.size;

        for(var i = 0; i < this.text.length; ++i) {
            var code = this.text.charCodeAt(i),
                ch = this.text.charAt(i);

            //if this is a newline
            if(/(?:\r\n|\r|\n)/.test(ch)) {
                lineWidths.push(pos.x);
                maxLineWidth = Math.max(maxLineWidth, pos.x);
                line++;

                pos.x = 0;
                pos.y += font.lineHeight;
                prevCode = null;
                continue;
            }

            //get this character's data
            var data = font.chars[code];

            if(!data) continue;

            //apply kernings
            if(prevCode && data[prevCode]) {
                pos.x += data.kerning[prevCode] || 0;
            }

            //add character
            chars.push({
                texture: data.texture,
                line: line,
                code: code,
                x: pos.x + data.xOffset,
                y: pos.y + data.yOffset
            });

            //advance the position
            pos.x += (this.monospace || data.xAdvance);

            //remember this code for kernings next char
            prevCode = code;
        }

        //final width
        lineWidths.push(pos.x);
        maxLineWidth = Math.max(maxLineWidth, pos.x);

        //unfortunately to do alignment, we have to loop through the lines to get
        //the offsets we need, then loop through characters to apply it. If we didn't
        //support alignment, then characters could be drawn in the above loop, but nooo...
        var lineAlignOffsets = [],
            align = this._style.align,
            offset = 0;

        for(i = 0; i <= line; ++i) {
            offset = 0;
            if(align === 'right')
                offset = maxLineWidth - lineWidths[i];
            else if(align === 'center')
                offset = (maxLineWidth - lineWidths[i]) / 2;

            lineAlignOffsets.push(offset);
        }

        //now add each character
        var lenChars = chars.length,
            lenChildren = this.children.length,
            tint = this._style.tint || 0xFFFFFF,
            child;

        for(i = 0; i < lenChars; ++i) {
            child = i < lenChildren ? this.children[i] : this.sprites.create(chars[i].texture);

            child.setTexture(chars[i].texture);

            child.position.x = (chars[i].x + lineAlignOffsets[chars[i].line]) * scale;
            child.position.y = chars[i].y * scale;
            child.scale.x = child.scale.y = scale;
            child.tint = tint;

            if(!child.parent)
                this.addChild(child);
        }

        //remove unnecesary children and free into pool
        while(this.children.length > lenChars) {
            child = this.children[this.children.length - 1];
            this.sprites.free(child);
            this.removeChild(child);
        }

        //set the width/height
        this.width = maxLineWidth * scale;
        this.height = (pos.y + font.lineHeight) * scale;
    },
    /**
     * Clones this BitmapText to get another just like it
     *
     * @method clone
     * @return BitmapText
     */
    clone: function() {
        return new BitmapText(this._text, this.font, this._style);
    },
    /**
     * Updates the text when dirty, called each frame by PIXI's render methods
     *
     * @method updateTransform
     * @private
     */
    updateTransform: function() {
        if(this.dirty) {
            this.renderText();

            this.dirty = false;
        }

        SpriteBatch.prototype.updateTransform.call(this);
    }
});

/**
 * The text that will be rendered.
 *
 * @property text
 * @type String
 */
Object.defineProperty(BitmapText.prototype, 'text', {
    get: function() {
        return this._text;
    },
    set: function(text) {
        this._text = text;
        this.dirty = true;
    }
});

/**
 * The style of the text to be rendered. Valid properties are
 * `size` and `align`.
 *
 * @property style
 * @type Object
 */
Object.defineProperty(BitmapText.prototype, 'style', {
    get: function() {
        return this._style;
    },
    set: function(style) {
        this._style.size = (style && style.size) || this._style.size;
        this._style.align = (style && style.align) || this._style.align;
        this.dirty = true;
    }
});

/**
 * The size of the text to render.
 *
 * @property size
 * @type Number
 */
Object.defineProperty(BitmapText.prototype, 'size', {
    get: function() {
        return this._style.size;
    },
    set: function(size) {
        this._style.size = size || this._style.size;
        this.dirty = true;
    }
});

/**
 * The alignment of the text to render, valid values are `'left'`,
 * `'right'`, or `'center'`.
 *
 * @property align
 * @type String
 */
Object.defineProperty(BitmapText.prototype, 'align', {
    get: function() {
        return this._style.align;
    },
    set: function(align) {
        this._style.align = align || this._style.align;
        this.dirty = true;
    }
});

/**
 * Parses an XML font file into a font object that can be passed into a BitmapText instance.
 * This is called by the Cache when XML bitmap data is added.
 *
 * @method parseXML
 * @param key {String} The cache key for the font
 * @param xml {Document} The XML document to parse
 * @param texture {Texture} The texture to use for creating character textures
 * @static
 */
BitmapText.parseXMLFont = function(key, xml, texture) {
    var btx = texture.baseTexture;

    if(!xml.getElementsByTagName('font')) {
        utils.warn('Invalid XML for BitmapText.parseXML(), missing <font> tag. Full XML:', xml);
    }

    var data = {},
        info = xml.getElementsByTagName('info')[0],
        common = xml.getElementsByTagName('common')[0];

    data.name = info.attributes.getNamedItem('face').nodeValue;
    data.size = parseInt(info.attributes.getNamedItem('size').nodeValue, 10);
    data.lineHeight = parseInt(common.attributes.getNamedItem('lineHeight').nodeValue, 10);
    data.chars = {};

    //parse characters
    var chars = xml.getElementsByTagName('char');

    for(var i = 0, il = chars.length; i < il; ++i) {
        var letter = chars[i],
            attrs = letter.attributes,
            code = parseInt(attrs.getNamedItem('id').nodeValue, 10),
            rect = new Rectangle(
                parseInt(attrs.getNamedItem('x').nodeValue, 10),
                parseInt(attrs.getNamedItem('y').nodeValue, 10),
                parseInt(attrs.getNamedItem('width').nodeValue, 10),
                parseInt(attrs.getNamedItem('height').nodeValue, 10)
            ),
            tx = PIXI.TextureCache[key + '_' + code] = new Texture(btx, rect);

        data.chars[code] = {
            xOffset: parseInt(attrs.getNamedItem('xoffset').nodeValue, 10),
            yOffset: parseInt(attrs.getNamedItem('yoffset').nodeValue, 10),
            xAdvance: parseInt(attrs.getNamedItem('xadvance').nodeValue, 10),
            kerning: {},
            texture: tx
        };
    }

    //parse kernings
    var kernings = xml.getElementsByTagName('kerning');

    for(i = 0, il = kernings.length; i < il; ++i) {
        var kern = kernings[i],
            attrs2 = kern.attributes,
            first = parseInt(attrs2.getNamedItem('first').nodeValue, 10),
            second = parseInt(attrs2.getNamedItem('second').nodeValue, 10),
            amount = parseInt(attrs2.getNamedItem('amount').nodeValue, 10);

        data.chars[second].kerning[first] = amount;
    }

    PIXI.BitmapText.fonts[data.name] = data;

    return data;
};

module.exports = BitmapText;

},{"../display/Sprite":19,"../display/SpriteBatch":20,"../display/Texture":21,"../geom/Rectangle":37,"../math/Vector":47,"../utils/ObjectPool":66,"../utils/inherit":69,"../utils/utils":71,"pixi.js":6}],55:[function(_dereq_,module,exports){
/**
 * A Text Object will create a line(s) of text to split a line you can use "\n"
 * see <a href="http://www.goodboydigital.com/pixijs/docs/classes/Text.html">PIXI.Text</a>
 * for more information.
 *
 * @class Text
 * @extends Sprite
 * @constructor
 * @param text {String} The copy that you would like the text to display
 * @param [style] {Object} The style parameters
 * @param [style.font] {String} default "bold 20pt Arial" The style and size of the font
 * @param [style.fill="black"] {Object} A canvas fillstyle that will be used on the text eg "red", "#00FF00"
 * @param [style.align="left"] {String} An alignment of the multiline text ("left", "center" or "right")
 * @param [style.stroke] {String} A canvas fillstyle that will be used on the text stroke eg "blue", "#FCFF00"
 * @param [style.strokeThickness=0] {Number} A number that represents the thickness of the stroke. Default is 0 (no stroke)
 * @param [style.wordWrap=false] {Boolean} Indicates if word wrap should be used
 * @param [style.wordWrapWidth=100] {Number} The width at which text will wrap
 */
var Text = _dereq_('pixi.js').Text;

module.exports = Text;

},{"pixi.js":6}],56:[function(_dereq_,module,exports){
var Container = _dereq_('../display/Container'),
    Vector = _dereq_('../math/Vector'),
    Polygon = _dereq_('../geom/Polygon'),
    Ellipse = _dereq_('../geom/Ellipse'),
    Rectangle = _dereq_('../geom/Rectangle'),
    utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    math = _dereq_('../math/math'),
    PIXI = _dereq_('pixi.js');

/**
 * Tiled object group is a special layer that contains entities
 * TODO: This is all trash
 *
 * @class ObjectGroup
 * @extends Container
 * @constructor
 * @param map {Tilemap} The tilemap instance that this belongs to
 * @param group {Object} All the settings for the layer
 */
var ObjectGroup = function(map, group) {
    Container.call(this, group);

    /**
     * The map instance this object group belongs to
     *
     * @property map
     * @type Tilemap
     */
    this.map = map;

    /**
     * The game instance this object group belongs to
     *
     * @property game
     * @type Game
     */
    this.game = map.game;

    /**
     * The state instance this object group belongs to
     *
     * @property state
     * @type Game
     */
    this.state = map.state;

    /**
     * The name of the group
     *
     * @property name
     * @type String
     * @default ''
     */
    this.name = group.name || '';

    /**
     * The color to display objects in this group
     *
     * @property color
     * @type
     */
    this.color = group.color;

    /**
     * The user-defined properties of this group. Usually defined in the TiledEditor
     *
     * @property properties
     * @type Object
     */
    this.properties = group.properties || {};

    /**
     * The objects in this group that can be spawned
     *
     * @property objects
     * @type Array
     */
    this.objects = group.objects;

    //translate some tiled properties to our inherited properties
    /**
     * The type of the layer, should always be 'objectgroup'
     *
     * @property type
     * @type String
     * @default 'objectgroup'
     */
    this.type = group.type || 'objectgroup';

    //translate some tiled properties to our inherited properties
    this.position.x = group.x || 0;
    this.position.y = group.y || 0;
    this.alpha = group.opacity !== undefined ? group.opacity : 1;
    this.visible = group.visible !== undefined ? group.visible : true;
};

inherit(ObjectGroup, Container, {
    getBounds: function() {
        return this.map.getBounds();
    },
    /**
     * Spawns all the entities associated with this layer, and properly sets their attributes
     *
     * @method spawn
     * @return {ObjectGroup} Returns itself.
     * @chainable
     */
    spawn: function() {
        var game = this.game; //this.Tilemap.GameState.Game

        //we go through these backwards so that things that are higher in the
        //list of object gets rendered on top.
        for(var i = this.objects.length - 1; i >= 0; --i) {
            var o = this.objects[i],
                props = utils.parseTiledProperties(o.properties) || {},
                set,
                interactive,
                obj;

            props.tileprops = {};

            //create a sprite with that texture
            if(o.gid) {
                set = this.parent.getTileset(o.gid);

                if(set) {
                    props.texture = set.getTileTexture(o.gid);
                    props.tileprops = set.getTileProperties(o.gid);

                    //if no hitArea then use the tileset's if available
                    if(!props.hitArea) {
                        if(props.tileprops.hitArea)
                            props.hitArea = props.tileprops.hitArea;
                        else
                            props.hitArea = set.properties.hitArea;
                    }
                }
            }
            //non-sprite object (usually to define an "area" on a map)
            else {
                if(!props.hitArea) {
                    //define a hitArea
                    if(o.polyline)
                        props.hitArea = this._getPolyline(o);
                    else if(o.polygon)
                        props.hitArea = this._getPolygon(o);
                    else if(o.ellipse)
                        props.hitArea = this._getEllipse(o);
                    else
                        props.hitArea = this._getRectangle(o);
                }
            }

            o.name = o.name || props.name || props.tileprops.name;
            o.type = o.type || props.type || props.tileprops.type;

            //a manually specified string texture
            if(typeof props.texture === 'string') {
                props.texture = game.cache.getTexture(props.texture);
            }

            //just a regular DisplayObject
            if(!props.texture) {
                obj = new Container();

                obj.width = o.width;
                obj.height = o.height;
                obj.name = o.name;
                obj.type = o.type;
                obj.rotation = o.rotation;

                //these are treated as sensor bodies, so always enable physics
                obj.position.x = o.x;
                obj.position.y = o.y;

                obj.sensor = true;
                obj.hitArea = props.hitArea;

                obj.enablePhysics(game.physics);
            } else {
                //some variable for the user if they want them
                //these will be passed through to a custom sprite if wanted
                props.width = o.width;
                props.height = o.height;
                props.zIndex = this.zIndex;

                obj = game.spritepool.create(o.name, props.texture, props);

                //assign some values
                obj.name = o.name;
                obj.type = o.type;
                obj.position.x = o.x;
                obj.position.y = o.y;

                obj.mass = props.mass || props.tileprops.mass;
                obj.inertia = props.inertia || props.tileprops.inertia;
                obj.friction = props.friction || props.tileprops.friction;
                obj.sensor = props.sensor || props.tileprops.sensor;
                obj.hitArea = props.hitArea;
                obj.blendMode = (props.blendMode || this.properties.blendMode) ? PIXI.blendModes[(props.blendMode || this.properties.blendMode)] : PIXI.blendModes.NORMAL;

                var a = props.anchor || props.tileprops.anchor;
                obj.anchor.y = a ? a[1] : 1;
                obj.anchor.x = a ? a[0] : (this.parent.orientation === 'isometric' ? 0.5 : 0);

                if(obj.mass) {
                    obj.enablePhysics(game.physics);
                }

                if(props.tileprops) {
                    if(props.tileprops.flippedX) {
                        obj.scale.x = -1;
                        obj.anchor.x = a ? a[0] : 1;
                    }

                    if(props.tileprops.flippedY) {
                        obj.scale.y = -1;
                        obj.anchor.y = a ? a[1] : 0;
                    }

                    //IDK if this is the correct angle, there are no docs for `rotatedCW`
                    if(props.tileprops.rotatedCW) {
                        obj.rotation = math.degreesToRadians(45);
                    }
                }

                if(props.animation || props.tileprops.animation) {
                    if(obj.goto) {
                        obj.goto(0, props.animation || props.tileprops.animation).play();
                    }
                }

                //set some more stuffz
                if(typeof o.rotation === 'number')
                    obj.setRotation(o.rotation);
            }

            //visible was recently added to Tiled, default old versions to true
            obj.visible = o.visible !== undefined ? !!o.visible : true;

            if(this.parent.orientation === 'isometric') {
                var toTileX = o.x / this.parent.tileSize.x,
                    toTileY = o.y / this.parent.tileSize.y;

                //This cannot be the simplest form of this...
                o.x = (toTileX * this.parent.tileSize.x) - ((toTileY - 1) * (this.parent.tileSize.x / 2));
                o.y = (toTileY * this.parent.tileSize.y / 2) + (toTileX * this.parent.tileSize.y);
            }

            interactive = this._getInteractive(set, props);

            //pass through all events
            if(interactive) {
                obj.interactive = interactive;

                obj.click = this.onObjectEvent.bind(this, 'click', obj);
                obj.mousedown = this.onObjectEvent.bind(this, 'mousedown', obj);
                obj.mouseup = this.onObjectEvent.bind(this, 'mouseup', obj);
                obj.mousemove = this.onObjectEvent.bind(this, 'mousemove', obj);
                obj.mouseout = this.onObjectEvent.bind(this, 'mouseout', obj);
                obj.mouseover = this.onObjectEvent.bind(this, 'mouseover', obj);
                obj.mouseupoutside = this.onObjectEvent.bind(this, 'mouseupoutside', obj);
            }

            //set custom properties
            obj.properties = {};
            for(var t in props.tileprops)
                obj.properties[t] = props.tileprops[t];
            for(var k in props)
                if(k !== 'tileprops')
                    obj.properties[k] = props[k];

            obj._objIndex = i;
            this.addChild(obj);
        }

        return this;
    },
    /**
     * Called internally whenever an event happens on an object, used to echo to the parent.
     *
     * @method onObjectEvent
     * @param eventName {String} The name of the event
     * @param obj {Container|Sprite} The object the event happened to
     * @param data {mixed} The event data that was passed along
     * @private
     */
    onObjectEvent: function(eventName, obj, data) {
        this.parent.onObjectEvent(eventName, obj, data);
    },
    /**
     * Creates a polygon from the vertices in a polygon Tiled property
     *
     * @method _getPolygon
     * @param obj {Object} The polygon Tiled object
     * @return {Polygon} The polygon created
     * @private
     */
    _getPolygon: function(o) {
        var points = [];
        for(var i = 0, il = o.polygon.length; i < il; ++i) {
            points.push(new Vector(o.polygon[i].x, o.polygon[i].y));
        }

        return new Polygon(points);
    },
    /**
     * Creates a polyline from the vertices in a polyline Tiled property
     *
     * @method _getPolyline
     * @param obj {Object} The polyline Tiled object
     * @return {Polygon} The polyline created
     * @private
     */
    _getPolyline: function(o) {
        var points = [];
        for(var i = 0, il = o.polyline.length; i < il; ++i) {
            points.push(new Vector(o.polyline[i].x, o.polyline[i].y));
        }

        return new Polygon(points);
    },
    /**
     * Creates a ellipse from the vertices in a ellipse Tiled property
     *
     * @method _getEllipse
     * @param obj {Object} The ellipse Tiled object
     * @return {Ellipse} The ellipse created
     * @private
     */
    _getEllipse: function(o) {
        return new Ellipse(0, 0, o.width, o.height);
    },
    /**
     * Creates a rectangle from the vertices in a rectangle Tiled property
     *
     * @method _getRectangle
     * @param obj {Object} The rectangle Tiled object
     * @return {Rectangle} The rectangle created
     * @private
     */
    _getRectangle: function(o) {
        return new Rectangle(0, 0, o.width, o.height);
    },
    /**
     * Checks if an object should be marked as interactive
     *
     * @method _getInteractive
     * @param set {Tileset} The tileset for the object
     * @param props {Object} The Tiled properties object
     * @return {Boolean} Whether or not the item is interactive
     * @private
     */
    _getInteractive: function(set, props) {
        //TODO: This is wrong, if 'false' is set on a lower level a higher level will override
        //first check the lowest level value (on the tile iteself)
        return props.interactive || //obj interactive
                props.tileprops.interactive || //tile object interactive
                (set && set.properties.interactive) || //tileset interactive
                this.properties.interactive || //layer interactive
                this.parent.properties.interactive; //map interactive
    },
    /**
     * Despawns all the sprites associated with this layer
     *
     * @method despawn
     * @return {ObjectGroup} Returns itself.
     * @chainable
     */
    despawn: function() {
        return Container.prototype.destroyAllChildren.call(this);
    },
    /**
     * Destroys the group completely
     *
     * @method destroy
     */
    destroy: function() {
        this.despawn();
        Container.prototype.destroy.call(this);

        this.map = null;
        this.game = null;
        this.state = null;
        this.name = null;
        this.color = null;
        this.properties = null;
        this.objects = null;
        this.type = null;
    }
});

module.exports = ObjectGroup;

},{"../display/Container":16,"../geom/Ellipse":35,"../geom/Polygon":36,"../geom/Rectangle":37,"../math/Vector":47,"../math/math":48,"../utils/inherit":69,"../utils/utils":71,"pixi.js":6}],57:[function(_dereq_,module,exports){
var Sprite = _dereq_('../display/Sprite'),
    inherit = _dereq_('../utils/inherit');

/**
 * Base Tile implementation, a tile is a single tile in a tilemap layer
 *
 * @class Tile
 * @extends Sprite
 * @constructor
 * @param texture {Texture} The texture of the tile
 */
var Tile = function(texture) {
    //call base ctor
    Sprite.call(this, texture);
};

inherit(Tile, Sprite, {
});

module.exports = Tile;

},{"../display/Sprite":19,"../utils/inherit":69}],58:[function(_dereq_,module,exports){
var SpriteBatch = _dereq_('../display/SpriteBatch'),
    Rectangle = _dereq_('../geom/Rectangle'),
    Vector = _dereq_('../math/Vector'),
    Texture = _dereq_('../display/Texture'),
    Tile = _dereq_('./Tile'),
    math = _dereq_('../math/math'),
    utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    support = _dereq_('../utils/support'),
    PIXI = _dereq_('pixi.js');

/**
 * The Tilelayer is the visual tiled layer that actually displays on the screen
 *
 * This class will be created by the Tilemap, there shouldn't be a reason to
 * create an instance on your own.
 *
 * @class Tilelayer
 * @extends SpriteBatch
 * @constructor
 * @param map {Tilemap} The tilemap instance that this belongs to
 * @param layer {Object} All the settings for the layer
 */
//see: https://github.com/GoodBoyDigital/pixi.js/issues/48
var Tilelayer = function(map, layer) {
    SpriteBatch.call(this);

    /**
     * The map instance this tilelayer belongs to
     *
     * @property map
     * @type Tilemap
     */
    this.map = map;

    /**
     * The state instance this tilelayer belongs to
     *
     * @property state
     * @type Game
     */
    this.state = map.state;

    /**
     * The state instance this tilelayer belongs to
     *
     * @property state
     * @type Game
     */
    this.state = map.state;

    /**
     * The current map of all tiles on the screen
     *
     * @property tiles
     * @type Object
     */
    this.tiles = [];

    /**
     * The name of the layer
     *
     * @property name
     * @type String
     * @default ''
     */
    this.name = layer.name || '';

    /**
     * The size of the layer
     *
     * @property size
     * @type Vector
     * @default new Vector(1, 1)
     */
    this.size = new Vector(layer.width || 0, layer.height || 0);

    /**
     * The tile IDs of the tilemap
     *
     * @property tileIds
     * @type Uint32Array
     */
    this.tileIds = support.typedArrays ? new Uint32Array(layer.data) : layer.data;

    /**
     * The user-defined properties of this group. Usually defined in the TiledEditor
     *
     * @property properties
     * @type Object
     */
    this.properties = utils.parseTiledProperties(layer.properties) || {};

    /**
     * The Tiled type of tile layer, should always be 'tilelayer'
     *
     * @property type
     * @type String
     * @default 'tilelayer'
     */
    this.type = layer.type || 'tilelayer';

    /**
     * Is this layer supposed to be preRendered?
     *
     * @property preRender
     * @type Boolean
     * @default false
     */
    this.preRender = layer.preRender || this.properties.preRender || this.map.properties.preRender || false;

    /**
     * The size of a chunk when pre rendering
     *
     * @property chunkSize
     * @type Vector
     * @default new Vector(512, 512)
     */
    this.chunkSize = new Vector(
        layer.chunkSizeX || layer.chunkSize || this.properties.chunkSizeX || this.properties.chunkSize || 512,
        layer.chunkSizeY || layer.chunkSize || this.properties.chunkSizeY || this.properties.chunkSize || 512
    );

    //translate some tiled properties to our inherited properties
    this.position.x = layer.x || 0;
    this.position.y = layer.y || 0;
    this.alpha = layer.opacity !== undefined ? layer.opacity : 1;
    this.visible = layer.visible !== undefined ? layer.visible : true;

    //some private trackers
    this._preRendered = false;
    this._tilePool = [];
    this._buffered = { left: false, right: false, top: false, bottom: false };
    this._panDelta = new Vector();
    this._rendered = new Rectangle();

    this.physicsContainer = new SpriteBatch();
    this.createPhysicalTiles();
};

inherit(Tilelayer, SpriteBatch, {
    getBounds: function() {
        return this.map.getBounds();
    },
    createPhysicalTiles: function() {
        var tid, tex, set, props, tile,
            szx = this.map.size.x,
            tsx = this.map.tileSize.x,
            tsy = this.map.tileSize.y;

        for(var i = 0; i < this.tileIds.length; ++i) {
            tid = this.tileIds[i];
            set = this.map.getTileset(tid);

            if(!set) continue;

            props = set.getTileProperties(tid);

            if(!props.mass) continue;

            tex = set.getTileTexture(tid);
            tile = new Tile(tex);
            this.physicsContainer.addChild(tile);

            tile.mass = props.mass;
            tile.hitArea = props.hitArea || set.properties.hitArea;
            tile.setPosition(
                ((i % szx) * tsx) + set.tileoffset.x,
                (math.floor(i / szx) * tsy) + set.tileoffset.y + tsy
            );

            tile.enablePhysics(this.state.physics);
        }
    },
    /**
     * Creates all the tile sprites needed to display the layer
     *
     * @method resize
     * @param width {Number} The number of tiles in the X direction to render
     * @param height {Number} The number of tiles in the Y direction to render
     * @return {Tilelayer} Returns itself.
     * @chainable
     */
    render: function(x, y, width, height) {
        if(this.preRender) {
            if(!this._preRendered) {
                this._preRender();
            } else {
                for(var c = this.children.length - 1; c > -1; --c) {
                    this.children[c].visible = true;
                }
            }

            return;
        }

        //copy down our tilesize
        if(!this.tileSize)
            this.tileSize = this.map.tileSize;

        //clear all the visual tiles
        this.clearTiles();

        //render the tiles on the screen
        this._renderTiles(x, y, width, height);

        return this;
    },
    /**
     * Renders the map onto different canvases, one per chunk. This only runs once
     * then the canvases are used as a textures for tiles the size of chunks.
     *
     * @method _preRender
     * @private
     */
    _preRender: function() {
        if(!this.visible)
            return;

        this._preRendered = true;
        this.tileSize = this.chunkSize.clone();

        var world = this.map,
            width = world.size.x * world.tileSize.x,
            height = world.size.y * world.tileSize.y,
            xChunks = math.ceil(width / this.chunkSize.x),
            yChunks = math.ceil(height / this.chunkSize.y);

        //for each chunk
        for(var x = 0; x < xChunks; ++x) {
            for(var y = 0; y < yChunks; ++y) {
                var cw = (x === xChunks - 1) ? width - (x * this.chunkSize.x) : this.chunkSize.x,
                    ch = (y === yChunks - 1) ? height - (y * this.chunkSize.y) : this.chunkSize.y;

                this._preRenderChunk(x, y, cw, ch);
            }
        }
    },
    /**
     * Renders a single chunk to a single canvas and creates/places the tile instance for it.
     *
     * @method _preRenderChunk
     * @param cx {Number} The x-coord of this chunk's top left
     * @param cy {Number} The y-coord of this chunk's top left
     * @param w {Number} The width of this chunk
     * @param h {Number} The height of this chunk
     * @private
     */
    _preRenderChunk: function(cx, cy, w, h) {
        var world = this.map,
            tsx = world.tileSize.x,
            tsy = world.tileSize.y,
            xTiles = w / tsx,
            yTiles = h / tsy,
            nx = (cx * this.chunkSize.x) % tsx,
            ny = (cy * this.chunkSize.y) % tsy,
            tx = math.floor(cx * this.chunkSize.x / tsx),
            ty = math.floor(cy * this.chunkSize.y / tsy),
            sx = world.size.x,
            sy = world.size.y,
            canvas = document.createElement('canvas'),
            ctx = canvas.getContext('2d');

        canvas.width = w;
        canvas.height = h;

        //draw all the tiles in this chunk to the canvas
        for(var x = 0; x < xTiles; ++x) {
            for(var y = 0; y < yTiles; ++y) {
                if(x + tx < sx && y + ty < sy) {
                    var id = ((x + tx) + ((y + ty) * sx)),
                        tid = this.tileIds[id],
                        set = world.getTileset(tid),
                        tex, frame;

                    if(set) {
                        tex = set.getTileTexture(tid);
                        frame = tex.frame;

                        ctx.drawImage(
                            tex.baseTexture.source,
                            frame.x,
                            frame.y,
                            frame.width,
                            frame.height,
                            (x * tsx) - nx + set.tileoffset.x,
                            (y * tsy) - ny + set.tileoffset.y,
                            frame.width,
                            frame.height
                        );
                    }
                }
            }
        }

        //use the canvas as a texture for a tile to display
        var tile = new Tile(Texture.fromCanvas(canvas));
        tile.setPosition(
            cx * this.chunkSize.x,
            cy * this.chunkSize.y
        );

        if(!this.tiles[cx])
            this.tiles[cx] = {};

        this.addChild(tile);
        this.tiles[cx][cy] = tile;
    },
    /**
     * Renders the tiles for the viewport
     *
     * @method _renderTiles
     * @param sx {Number} The x-coord in the map to start rendering
     * @param sy {Number} The y-coord in the map to start rendering
     * @param sw {Number} The width of the viewport
     * @param sh {Number} The height of the viewport
     * @private
     */
    _renderTiles: function(sx, sy, sw, sh) {
        //convert to tile coords
        sx = math.floor(sx / this.map.scaledTileSize.x);
        sy = math.floor(sy / this.map.scaledTileSize.y);

        //ensure we don't go below 0
        sx = sx < 0 ? 0 : sx;
        sy = sy < 0 ? 0 : sy;

        //convert to tile sizes
        sw = math.ceil(sw / this.map.scaledTileSize.x) + 1;
        sh = math.ceil(sh / this.map.scaledTileSize.y) + 1;

        //ensure we don't go outside the map size
        sw = (sx + sw > this.map.size.x) ? (this.map.size.x - sx) : sw;
        sh = (sy + sh > this.map.size.y) ? (this.map.size.y - sy) : sh;

        //render new sprites
        var endX = sx + sw,
            endY = sy + sh;

        for(var x = sx; x < endX; ++x) {
            for(var y = sy; y < endY; ++y) {
                this.moveTileSprite(-1, -1, x, y);
            }
        }

        //set rendered area
        this._rendered.x = sx;
        this._rendered.y = sy;
        this._rendered.width = sw - 1;
        this._rendered.height = sh - 1;

        //reset buffered status
        this._buffered.left = this._buffered.right = this._buffered.top = this._buffered.bottom = false;

        //reset panDelta
        this._panDelta.x = this.state.world.position.x % this.map.scaledTileSize.x;
        this._panDelta.y = this.state.world.position.y % this.map.scaledTileSize.y;
    },
    /**
     * Frees a tile in the list back into the pool
     *
     * @method _freeTile
     * @param tx {Number} The x-coord of the tile in tile coords (not world coords)
     * @param ty {Number} The y-coord of the tile in tile coords (not world coords)
     * @private
     */
    _freeTile: function(tx, ty) {
        if(this.tiles[tx] && this.tiles[tx][ty]) {
            this.clearTile(this.tiles[tx][ty]);
            this.tiles[tx][ty] = null;
        }
    },
    /**
     * Clears all the tiles currently used to render the layer
     *
     * @method clearTiles
     * @param remove {Boolean} Should this tile be completely removed (never to bee seen again)
     * @return {Tilelayer} Returns itself.
     * @chainable
     */
    clearTiles: function(remove) {
        var c;

        if(this.preRender && !remove) {
            for(c = this.children.length - 1; c > -1; --c) {
                this.children[c].visible = false;
            }

            return;
        }

        //force rerender later
        this._preRendered = false;

        for(c = this.children.length - 1; c > -1; --c) {
            this.clearTile(this.children[c], remove);
        }

        this.tiles.length = 0;

        return this;
    },
    /**
     * Clears a tile currently used to render the layer
     *
     * @method clearTile
     * @param tile {Tile} The tile object to clear
     * @param remove {Boolean} Should this tile be completely removed (never to bee seen again)
     * @return {Tilelayer} Returns itself.
     * @chainable
     */
    clearTile: function(tile, remove) {
        tile.visible = false;
        //tile.disablePhysics();

        if(remove)
            this.removeChild(tile);
        else
            this._tilePool.push(tile);

        return this;
    },
    /**
     * Moves a tile sprite from one position to another, and creates a new tile
     * if the old position didn't have a sprite
     *
     * @method moveTileSprite
     * @param fromTileX {Number} The x coord of the tile in units of tiles (not pixels) to move from
     * @param fromTileY {Number} The y coord of the tile in units of tiles (not pixels) to move from
     * @param toTileX {Number} The x coord of the tile in units of tiles (not pixels) to move to
     * @param toTileY {Number} The y coord of the tile in units of tiles (not pixels) to move to
     * @return {Tile} The sprite to display
     */
    moveTileSprite: function(fromTileX, fromTileY, toTileX, toTileY) {
        //free the tiles we are dealing with
        this._freeTile(toTileX, toTileY);
        this._freeTile(fromTileX, fromTileY);

        //if off the map, just ignore it
        if(toTileX < 0 || toTileY < 0 || toTileX >= this.map.size.x || toTileY >= this.map.size.y) {
            return;
        }

        var tile,
            id = (toTileX + (toTileY * this.size.x)),
            tileId = this.tileIds[id],
            set = this.map.getTileset(tileId),
            texture,
            props,
            position,
            hitArea,
            interactive;

        //if no tileset, return
        if(!set) return;

        //grab some values for the tile
        texture = set.getTileTexture(tileId);
        props = set.getTileProperties(tileId);
        hitArea = props.hitArea || set.properties.hitArea;
        interactive = this._getInteractive(set, props);
        position = [
            (toTileX * this.map.tileSize.x) + set.tileoffset.x,
            (toTileY * this.map.tileSize.y) + set.tileoffset.y
        ];

        //due to the fact that we use top-left anchors for everything, but tiled uses bottom-left
        //we need to move the position of each tile down by a single map-tile height. That is why
        //there is an addition of "this.map.tileSize.y" to the coords
        position[1] += this.map.tileSize.y;

        //grab a new tile from the pool
        tile = this._tilePool.pop();

        //if we couldn't find a tile from the pool, then create a new tile
        if(!tile) {
            tile = new Tile(texture);
            tile.anchor.y = 1;
            this.addChild(tile);
        }

        //tile.collisionType = props.type;
        tile.interactive = interactive;
        tile.hitArea = hitArea;
        //tile.mass = props.mass || 0;
        tile.blendMode = (props.blendMode || this.properties.blendMode) ? PIXI.blendModes[(props.blendMode || this.properties.blendMode)] : PIXI.blendModes.NORMAL;

        tile.setTexture(texture);
        tile.setPosition(position[0], position[1]);
        tile.show();

        /*if(tile.mass) {
            tile.enablePhysics(this.state.physics);
        }*/

        //pass through all events
        if(interactive) {
            tile.click = this.onTileEvent.bind(this, 'click', tile);
            tile.mousedown = this.onTileEvent.bind(this, 'mousedown', tile);
            tile.mouseup = this.onTileEvent.bind(this, 'mouseup', tile);
            tile.mousemove = this.onTileEvent.bind(this, 'mousemove', tile);
            tile.mouseout = this.onTileEvent.bind(this, 'mouseout', tile);
            tile.mouseover = this.onTileEvent.bind(this, 'mouseover', tile);
            tile.mouseupoutside = this.onTileEvent.bind(this, 'mouseupoutside', tile);
        }

        //update sprite position in the map
        if(!this.tiles[toTileX])
            this.tiles[toTileX] = [];

        this.tiles[toTileX][toTileY] = tile;

        return tile;
    },
    /**
     * Called whenever a tile event occurs, this is used to echo to the parent.
     *
     * @method onTileEvent
     * @param eventName {String} The name of the event
     * @param tile {Tile} The tile the event happened to
     * @param data {mixed} The event data that was passed along
     * @private
     */
    onTileEvent: function(eventName, tile, data) {
        this.map.onTileEvent(eventName, tile, data);
    },
    /**
     * Checks if an object should be marked as interactive
     *
     * @method _getInteractive
     * @param set {Tileset} The tileset for the object
     * @param props {Object} The Tiled properties object
     * @return {Boolean} Whether or not the item is interactive
     * @private
     */
    _getInteractive: function(set, props) {
        //first check the lowest level value (on the tile iteself)
        return props.interactive || //obj interactive
                (set && set.properties.interactive) || //tileset interactive
                this.properties.interactive || //layer interactive
                this.map.properties.interactive; //map interactive
    },
    /**
     * Pans the layer around, rendering stuff if necessary
     *
     * @method pan
     * @param dx {Number|Point} The x amount to pan, if a Point is passed the dy param is ignored
     * @param dy {Number} The y ammount to pan
     * @return {Tilelayer} Returns itself.
     * @chainable
     */
    pan: function(dx, dy) {
        if(this.preRender)
            return;

        //track panning delta so we know when to render
        this._panDelta.x += dx;
        this._panDelta.y += dy;

        var tszX = this.map.scaledTileSize.x,
            tszY = this.map.scaledTileSize.y;

        //check if we need to build a buffer around the viewport
        //usually this happens on the first pan after a full render
        //caused by a viewport resize. We do this buffering here instead
        //of in the initial render because in the initial render, the buffer
        //may try to go negative which has no tiles. Plus doing it here
        //reduces the number of tiles that need to be created initially.

        //moving world right, so left will be exposed
        if(dx > 0 && !this._buffered.left)
            this._renderLeft(this._buffered.left = true);
        //moving world left, so right will be exposed
        else if(dx < 0 && !this._buffered.right)
            this._renderRight(this._buffered.right = true);

        //moving world down, so top will be exposed
        if(dy > 0 && !this._buffered.top)
            this._renderUp(this._buffered.top = true);
        //moving world up, so bottom will be exposed
        else if(dy < 0 && !this._buffered.bottom)
            this._renderDown(this._buffered.bottom = true);

        //Here is where the actual panning gets done, we check if the pan
        //delta is greater than a scaled tile and if so pan that direction.
        //The reason we do it in a while loop is because the delta can be
        //large than 1 scaled tile and may require multiple render pans
        //(this can happen if you can .pan(x, y) with large values)

        //moved position right, so render left
        while(this._panDelta.x >= tszX) {
            this._renderLeft();
            this._panDelta.x -= tszX;
        }

        //moved position left, so render right
        while(this._panDelta.x <= -tszX) {
            this._renderRight();
            this._panDelta.x += tszX;
        }

        //moved position down, so render up
        while(this._panDelta.y >= tszY) {
            this._renderUp();
            this._panDelta.y -= tszY;
        }

        //moved position up, so render down
        while(this._panDelta.y <= -tszY) {
            this._renderDown();
            this._panDelta.y += tszY;
        }
    },
    /**
     * Renders tiles to the left, pulling from the far right
     *
     * @method _renderLeft
     * @param [forceNew=false] {Boolean} If set to true, new tiles are created instead of trying to recycle
     * @private
     */
    _renderLeft: function(forceNew) {
        //move all the far right tiles to the left side
        for(var i = 0; i < this._rendered.height + 1; ++i) {
            this.moveTileSprite(
                forceNew ? -1 : this._rendered.right,
                forceNew ? -1 : this._rendered.top + i,
                this._rendered.left - 1,
                this._rendered.top + i
            );
        }
        this._rendered.x--;
        if(forceNew) this._rendered.width++;
    },
    /**
     * Renders tiles to the right, pulling from the far left
     *
     * @method _renderRight
     * @param [forceNew=false] {Boolean} If set to true, new tiles are created instead of trying to recycle
     * @private
     */
    _renderRight: function(forceNew) {
        //move all the far left tiles to the right side
        for(var i = 0; i < this._rendered.height + 1; ++i) {
            this.moveTileSprite(
                forceNew ? -1 : this._rendered.left,
                forceNew ? -1 : this._rendered.top + i,
                this._rendered.right + 1,
                this._rendered.top + i
            );
        }
        if(!forceNew) this._rendered.x++;
        if(forceNew) this._rendered.width++;
    },
    /**
     * Renders tiles to the top, pulling from the far bottom
     *
     * @method _renderUp
     * @param [forceNew=false] {Boolean} If set to true, new tiles are created instead of trying to recycle
     * @private
     */
    _renderUp: function(forceNew) {
        //move all the far bottom tiles to the top side
        for(var i = 0; i < this._rendered.width + 1; ++i) {
            this.moveTileSprite(
                forceNew ? -1 : this._rendered.left + i,
                forceNew ? -1 : this._rendered.bottom,
                this._rendered.left + i,
                this._rendered.top - 1
            );
        }
        this._rendered.y--;
        if(forceNew) this._rendered.height++;
    },
    /**
     * Renders tiles to the bottom, pulling from the far top
     *
     * @method _renderDown
     * @param [forceNew=false] {Boolean} If set to true, new tiles are created instead of trying to recycle
     * @private
     */
    _renderDown: function(forceNew) {
        //move all the far top tiles to the bottom side
        for(var i = 0; i < this._rendered.width + 1; ++i) {
            this.moveTileSprite(
                forceNew ? -1 : this._rendered.left + i,
                forceNew ? -1 : this._rendered.top,
                this._rendered.left + i,
                this._rendered.bottom + 1
            );
        }
        if(!forceNew) this._rendered.y++;
        if(forceNew) this._rendered.height++;
    },
    /**
     * Destroys the tile layer completely
     *
     * @method destroy
     */
    destroy: function() {
        SpriteBatch.prototype.destroy.call(this);

        this.clearTiles(true);

        this.state = null;
        this.name = null;
        this.size = null;
        this.tileIds = null;
        this.properties = null;
        this.type = null;
        this.position.x = null;
        this.position.y = null;
        this.alpha = null;
        this.visible = null;
        this.preRender = null;
        this.chunkSize = null;

        this._preRendered = null;
        this._tilePool = null;
        this._buffered = null;
        this._panDelta = null;
        this._rendered = null;
    }
});

module.exports = Tilelayer;

},{"../display/SpriteBatch":20,"../display/Texture":21,"../geom/Rectangle":37,"../math/Vector":47,"../math/math":48,"../utils/inherit":69,"../utils/support":70,"../utils/utils":71,"./Tile":57,"pixi.js":6}],59:[function(_dereq_,module,exports){
var Container = _dereq_('../display/Container'),
    ObjectGroup = _dereq_('./ObjectGroup'),
    Sprite = _dereq_('../display/Sprite'),
    Vector = _dereq_('../math/Vector'),
    Rectangle = _dereq_('../geom/Rectangle'),
    Tilelayer = _dereq_('./Tilelayer'),
    Tileset = _dereq_('./Tileset'),
    utils = _dereq_('../utils/utils'),
    math = _dereq_('../math/math'),
    inherit = _dereq_('../utils/inherit');

/**
 * Tiled map that represents an entire tile map with multiple layers or object groups.
 * Often it is easier to create a tilemap using the object factor on a world, rather
 * than doing it manually yourself.
 *
 * @class Tilemap
 * @extends Container
 * @constructor
 * @param state {State} The game state the map belongs to
 * @param map {Object} All the settings for the map
 * @param tilesetTextures {Object} An object whose keys are the tileset name,
 *      and whose values are the textures for the tileset. For example:
 *      `{ tileset1: new Texture(), ... }`
 */
var Tilemap = function(state, map, tilesetTextures) {
    //call base ctor
    Container.call(this, map);

    /**
     * The state instance this tilemap belongs to
     *
     * @property state
     * @type Game
     */
    this.state = state;

    /**
     * The game instance this tilemap belongs to
     *
     * @property game
     * @type Game
     */
    this.game = state.game;


    //Tiled Editor properties

    /**
     * The user-defined properties
     *
     * @property properties
     * @type Object
     * @default {}
     */
    this.properties = utils.parseTiledProperties(map.properties) || {};
    this.scale.x = this.properties.scale || 1;
    this.scale.y = this.properties.scale || 1;

    /**
     * The tile size
     *
     * @property tileSize
     * @type Vector
     */
    this.tileSize = new Vector(map.tilewidth, map.tileheight);

    /**
     * The size of the map
     *
     * @property size
     * @type Vector
     * @default new Vector(0, 0)
     */
    this.size = new Vector(map.width, map.height);

    /**
     * The orientation of the map
     *
     * @property orientation
     * @type String
     */
    this.orientation = map.orientation;

    /**
     * The version of the TMX format
     *
     * @property version
     * @type Number
     */
    this.version = map.version;

    /**
     * The background color of the map (since Tiled 0.9.0)
     *
     * @property backgroundColor
     * @type Number
     */
    this.backgroundColor = map.backgroundColor;

    //Custom/Optional properties

    /**
     * The tilesets used by this map
     *
     * @property tilesets
     * @type Array
     */
    this.tilesets = [];

    /**
     * The scaled tile size
     *
     * @property scaledTileSize
     * @type Vector
     */
    this.scaledTileSize = new Vector(
        map.tilewidth * this.scale.x,
        map.tileheight * this.scale.y
    );

    /**
     * The real size (size * scaledTileSize)
     *
     * @property realSize
     * @type Vector
     */
    this.realSize = new Vector(
        this.size.x * this.scaledTileSize.x,
        this.size.y * this.scaledTileSize.y
    );

    //create each tileset
    for(var t = 0, tl = map.tilesets.length; t < tl; ++t) {
        var ts = map.tilesets[t];
        this.tilesets.push(new Tileset(tilesetTextures[ts.name], ts));
    }

    //create each layer
    for(var i = 0, il = map.layers.length; i < il; ++i) {
        var lyr;

        switch(map.layers[i].type) {
        case 'tilelayer':
            lyr = new Tilelayer(this, map.layers[i]);
            break;

        case 'objectgroup':
            lyr = new ObjectGroup(this, map.layers[i]);
            break;

        case 'imagelayer':
            lyr = new Sprite(map.layers[i]);
            break;
        }

        this.addChild(lyr);
    }

    this._bounds = new Rectangle(0, 0, this.realSize.x, this.realSize.y);

    //update the world bounds
    var w = this.game.state.active.world;
    w.bounds.width = Math.max(w.bounds.width, this.realSize.x);
    w.bounds.height = Math.max(w.bounds.height, this.realSize.y);
};

inherit(Tilemap, Container, {
    getBounds: function() {
        this.scaledTileSize.set(
            this.tileSize.x * this.scale.x,
            this.tileSize.y * this.scale.y
        );
        this.realSize.set(
            this.size.x * this.scaledTileSize.x,
            this.size.y * this.scaledTileSize.y
        );
        this._bounds.width = math.min(this.realSize.x, this.state.game.width);
        this._bounds.height = math.min(this.realSize.y, this.state.game.height);

        return this._bounds;
    },
    /**
     * Gets the tileset that an ID is associated with
     *
     * @method getTileset
     * @param tileId {Number} The id of the tile to find the tileset for
     * @return {TiledTileset} Returns the tileset if found, undefined if not
     */
    getTileset: function(tileId) {
        for(var i = 0, il = this.tilesets.length; i < il; ++i)
            if(this.tilesets[i].contains(tileId))
                return this.tilesets[i];
    },
    /**
     * Destroys the tilemap instance
     *
     * @method destroy
     */
    destroy: function() {
        Container.prototype.destroy.call(this);

        this.game = null;
        this.properties = null;
        this.tileSize = null;
        this.size = null;
        this.orientation = null;
        this.version = null;
        this.backgroundColor = null;
        this.tilesets = null;
        this.scaledTileSize = null;
        this.realSize = null;
    },
    /**
     * Spawns all the objects in the ObjectGroups of this map
     *
     * @method spawnObjects
     * @return {Tilemap} Returns itself.
     * @chainable
     */
    spawnObjects: function() {
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.type === 'objectgroup') {
                o.spawn();
            }
        }

        return this;
    },
    /**
     * Spawns all the objects in the ObjectGroups of this map
     *
     * @method despawnObjects
     * @return {Tilemap} Returns itself.
     * @chainable
     */
    despawnObjects: function() {
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.type === 'objectgroup') {
                o.despawn();
            }
        }

        return this;
    },
    /**
     * Clears all the tiles that are currently used on all tile layers
     *
     * @method clearTiles
     * @return {Tilemap} Returns itself.
     * @chainable
     */
    clearTiles: function() {
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.type === 'tilelayer') {
                o.clearTiles();
            }
        }

        return this;
    },
    /**
     * Called by a Tilelayer when a tile event occurs. This is so you can listen for
     * the emitted events on the world instead of the tile itself.
     *
     * @method onTileEvent
     * @param eventName {String} The event name to emit, the prefix 'tile.' will be added to it
     * @param tile {Tile} The tile that has the event
     * @param data {InteractionData} The raw interaction object for the event
     * @private
     */
    onTileEvent: function(eventName, tile, data) {
        this.emit('tile.' + eventName, {
            tile: tile,
            data: data
        });
    },
    /**
     * Called by a ObjectGroup when an object event occurs. This is so you can listen for
     * the emitted events on the world instead of the tile itself.
     *
     * @method onObjectEvent
     * @param eventName {String} The event name to emit, the prefix 'object.' will be added to it
     * @param obj {Sprite|Container} The object that has the event
     * @param data {InteractionData} The raw interaction object for the event
     * @private
     */
    onObjectEvent: function(eventName, obj, data) {
        this.emit('object.' + eventName, {
            object: obj,
            data: data
        });
    },
    /**
     * Finds a layer based on the string name
     *
     * @method findLayer
     * @param name {String} The name of the layer to find
     * @return {Tilelayer|ObjectGroup|Sprite} Returns the layer if found, undefined if not
     */
    findLayer: function(name) {
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.name === name)
                return o;
        }
    },
    /**
     * Pans the map around
     *
     * @method pan
     * @param x {Number|Point} The x amount to pan, if a Point is passed the y param is ignored
     * @param y {Number} The y ammount to pan
     * @return {Tilemap} Returns itself.
     * @chainable
     */
    pan: function(x, y) {
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.pan)
                o.pan(x, y);
        }

        return this;
    },
    /**
     * Called on resize to render the viewport again
     *
     * @method render
     * @param x {Number} The x offset to consider the top-left
     * @param y {Number} The y offset to consider the top-left
     * @param width {Number} The width (in pixels) to render
     * @param height {Number} The height (in pixels) to render
     * @return {Tilemap} Returns itself.
     * @chainable
     */
    render: function(x, y, width, height) {
        //defaults
        x = x || -this.state.world.position.x;
        y = y || -this.state.world.position.y;
        width = width || this.game.width;
        height = height || this.game.height;

        //render the layers
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var o = this.children[i];

            if(o.render)
                o.render(x, y, width, height);
        }

        return this;
    }
});

Tilemap.parseXMLMap = function(data) {
    var mapElement = data.getElementsByTagName('map')[0],
        map = {
            version: parseFloat(mapElement.attributes.getNamedItem('version').nodeValue, 10),
            width: parseInt(mapElement.attributes.getNamedItem('width').nodeValue, 10),
            height: parseInt(mapElement.attributes.getNamedItem('height').nodeValue, 10),
            tilewidth: parseInt(mapElement.attributes.getNamedItem('tilewidth').nodeValue, 10),
            tileheight: parseInt(mapElement.attributes.getNamedItem('tileheight').nodeValue, 10),
            orientation: mapElement.attributes.getNamedItem('orientation').nodeValue,
            layers: [],
            properties: {},
            tilesets: []
        },
        i = 0,
        il = 0;

    //add the properties
    var mapprops = mapElement.getElementsByTagName('properties');
    for(i = 0, il = mapprops.length; i < il; ++i) {
        if(mapprops[i].parentNode === mapElement) {
            mapprops = mapprops.getElementsByTagName('property');

            for(var mp = 0; mp < mapprops.length; ++mp) {
                map.properties[mapprops[mp].attributes.getNamedItem('name').nodeValue] = mapprops[mp].attributes.getNamedItem('value').nodeValue;
            }

            break;
        }
    }

    //add the layers
    var layers = mapElement.childNodes;//getElementsByTagName('layer');

    for(i = 0, il = layers.length; i < il; ++i) {
        var node = layers[i];

        if(node.nodeName === 'layer') {
            var lyr = node,
                layer = {
                    type: 'tilelayer',
                    name: lyr.attributes.getNamedItem('name').nodeValue,
                    width: parseInt(lyr.attributes.getNamedItem('width').nodeValue, 10) || map.width,
                    height: parseInt(lyr.attributes.getNamedItem('height').nodeValue, 10) || map.height,
                    visible: lyr.attributes.getNamedItem('visible') ? lyr.attributes.getNamedItem('visible').nodeValue === '1' : true,
                    opacity: lyr.attributes.getNamedItem('opacity') ? parseFloat(lyr.attributes.getNamedItem('opacity').nodeValue, 10) : 1,
                    encoding: 'base64',
                    compression: '',
                    rawData: '',
                    data: '',
                    x: 0,
                    y: 0
                };

            //set encoding
            var dataElement = lyr.getElementsByTagName('data')[0];
            layer.encoding = dataElement.attributes.getNamedItem('encoding').nodeValue;

            //set data from the text node of the element
            layer.rawData = dataElement.firstChild.nodeValue.trim();

            //set compression
            if(dataElement.attributes.getNamedItem('compression')) {
                layer.compression = dataElement.attributes.getNamedItem('compression').nodeValue;

                var decomp = Tilemap._decompressBase64LayerData(layer.rawData, layer.encoding, layer.compression);

                layer.data = new Uint32Array(decomp.buffer, 0, decomp.length / 4);
            }
            //not compressed, just decode base64
            else if(layer.encoding === 'base64') {
                var decoded = Tilemap._decodeBase64LayerData(layer.rawData);

                layer.data = new Uint32Array(decoded.buffer, 0, decoded.length / 4);
            }

            map.layers.push(layer);
        } else if(node.nodeName === 'objectgroup') {
            var grp = node,
                group = {
                    type: 'objectgroup',
                    draworder: 'topdown',
                    name: grp.attributes.getNamedItem('name').nodeValue,
                    width: 0,
                    height: 0,
                    objects: [],
                    visible: grp.attributes.getNamedItem('visible') ? grp.attributes.getNamedItem('visible').nodeValue === '0' : true,
                    opacity: grp.attributes.getNamedItem('opacity') ? parseFloat(grp.attributes.getNamedItem('opacity').nodeValue, 10) : 1,
                    x: 0,
                    y: 0
                };

            var objects = grp.getElementsByTagName('object');
            for(var oj = 0; oj < objects.length; ++oj) {
                var obj = objects[oj],
                    object = {
                        gid: obj.attributes.getNamedItem('gid') ? parseInt(obj.attributes.getNamedItem('gid').nodeValue, 10) : null,
                        name: obj.attributes.getNamedItem('name') ? obj.attributes.getNamedItem('name').nodeValue : '',
                        type: obj.attributes.getNamedItem('type') ? obj.attributes.getNamedItem('type').nodeValue : '',
                        width: obj.attributes.getNamedItem('width') ? parseFloat(obj.attributes.getNamedItem('width').nodeValue, 10) : 0,
                        height: obj.attributes.getNamedItem('height') ? parseFloat(obj.attributes.getNamedItem('height').nodeValue, 10) : 0,
                        rotation: obj.attributes.getNamedItem('rotation') ? parseFloat(obj.attributes.getNamedItem('rotation').nodeValue, 10) : 0,
                        visible: obj.attributes.getNamedItem('visible') ? obj.attributes.getNamedItem('visible').nodeValue === '1' : true,
                        x: parseFloat(obj.attributes.getNamedItem('x').nodeValue, 10),
                        y: parseFloat(obj.attributes.getNamedItem('y').nodeValue, 10),
                        properties: {}
                    };

                if(object.gid === null)
                    delete object.gid;

                var props = obj.getElementsByTagName('properties');
                if(props.length) {
                    props = props.getElementsByTagName('property');
                    for(var pr = 0; pr < props.length; ++pr) {
                        object.properties[props[pr].attributes.getNamedItem('name').nodeValue] = props[pr].attributes.getNamedItem('value').nodeValue;
                    }
                }

                group.objects.push(object);
            }

            map.layers.push(group);
        }
    }

    //add the tilesets
    var tilesets = mapElement.getElementsByTagName('tileset');

    for(i = 0, il = tilesets.length; i < il; ++i) {
        var tset = tilesets[i],
            tiles = tset.getElementsByTagName('tile'),
            tileset = {
                name: tset.attributes.getNamedItem('name').nodeValue,
                firstgid: parseInt(tset.attributes.getNamedItem('firstgid').nodeValue, 10),
                tilewidth: parseInt(tset.attributes.getNamedItem('tilewidth').nodeValue, 10),
                tileheight: parseInt(tset.attributes.getNamedItem('tileheight').nodeValue, 10),
                margin: 0,
                spacing: 0,
                tileoffset: { x: 0, y: 0 },
                terrains: [],
                properties: {},
                tileproperties: {},
                tiles: {}
            };

        //add spacing / margin attributes if exist
        var spacing = tset.attributes.getNamedItem('spacing');
        if(spacing) {
            tileset.spacing = parseInt(spacing.nodeValue, 10);
        }

        var margin = tset.attributes.getNamedItem('margin');
        if(margin) {
            tileset.margin = parseInt(margin.nodeValue, 10);
        }

        //add .properties if element exists
        var tsetprops = tset.getElementsByTagName('properties');
        for(var tsp = 0; tsp < tsetprops.length; ++tsp) {
            if(tsetprops[tsp].parentNode === tset) {
                tsetprops = tsetprops.getElementsByTagName('property');

                if(tsetprops.length) {
                    for(var p = 0; p < tsetprops.length; ++p) {
                        var tsetprop = tsetprops[p];

                        tileset.properties[tsetprop.attributes.getNamedItem('name').nodeValue] = tsetprop.attributes.getNamedItem('value').nodeValue;
                    }
                }

                break;
            }
        }

        //add .tiles if multi-image set
        for(var t = 0; t < tiles.length; ++t) {
            var tile = tiles[t],
                id = tile.attributes.getNamedItem('id').nodeValue,
                img = tile.getElementsByTagName('image');

            tileset.tiles[id] = {};

            //add attributes into the object
            for(var ta = 0; ta < tile.attributes.length; ++ta) {
                var tileatr = tile.attributes[ta];
                if(tileatr.name === 'id') continue;

                switch(tileatr.name) {
                    case 'terrain':
                        tileset.tiles[id].terrain = tileatr.value.sply(',');
                        break;

                    case 'probability':
                        tileset.tiles[id].probability = parseFloat(tileatr.value, 10);
                        break;
                }
            }

            //check if it has an image child
            if(img.length) {
                tileset.tiles[id] = tileset.tiles[id] || {};
                tileset.tiles[id].image = img[0].attributes.getNamedItem('source').nodeValue;
            }

            //add all the tile properties
            var tileprops = tile.getElementsByTagName('properties');
            if(tileprops.length) {
                tileset.tileproperties[id] = {};
                tileprops = tileprops[0].getElementsByTagName('property');
                for(var tp = 0; tp < tileprops.length; ++tp) {
                    var tileprop = tileprops[tp];
                    tileset.tileproperties[id][tileprop.attributes.getNamedItem('name').nodeValue] = tileprop.attributes.getNamedItem('value').nodeValue;
                }
            }
        }

        //check for terraintypes and add those
        var terrains = tset.getElementsByTagName('terraintypes');
        if(terrains.length) {
            terrains = terrains[0].getElementsByTagName('terrain');
            for(var tr = 0; tr < terrains.length; ++tr) {
                tileset.terrains.push({
                    name: terrains[tr].attributes.getNamedItem('name').nodeValue,
                    tile: parseInt(terrains[tr].attributes.getNamedItem('tile').nodeValue, 10)
                });
            }
        }

        //check for tileoffset and add that
        var offset = tset.getElementsByTagName('tileoffset');
        if(offset.length) {
            tileset.tileoffset.x = parseInt(offset[0].attributes.getNamedItem('x').nodeValue, 10);
            tileset.tileoffset.y = parseInt(offset[0].attributes.getNamedItem('y').nodeValue, 10);
        }

        //add image, imagewidth, imageheight
        var image = tset.getElementsByTagName('image');
        if(image.length === 1 && image[0].parentNode === tset) {
            tileset.image = image[0].attributes.getNamedItem('source').nodeValue;
            tileset.imagewidth = parseInt(image[0].attributes.getNamedItem('width').nodeValue, 10);
            tileset.imageheight = parseInt(image[0].attributes.getNamedItem('height').nodeValue, 10);
        }

        map.tilesets.push(tileset);
    }

    return map;
};

Tilemap._decodeBase64LayerData = function(raw) {
    return new Uint8Array(window.atob(raw).split('').map(Tilemap._mapCharactersToCodes));
};

Tilemap._mapCharactersToCodes = function(ch) {
    return ch.charCodeAt(0);
};

Tilemap._decompressBase64LayerData = function(raw, encoding, compression) {
    //TODO: This assumes base64 encoding
    var zlib = _dereq_('zlibjs'),
        data = Tilemap._decodeBase64LayerData(raw),
        decomp;

    //decompress
    switch(compression) {
        case 'gzip':
            decomp = new zlib.gunzipSync(data);
            break;

        case 'zlib':
            decomp = new zlib.inflateSync(data);
            break;
    }

    return decomp;
};

/*
Tilemap.fromCSV = function(state, data) {

};
*/

module.exports = Tilemap;

},{"../display/Container":16,"../display/Sprite":19,"../geom/Rectangle":37,"../math/Vector":47,"../math/math":48,"../utils/inherit":69,"../utils/utils":71,"./ObjectGroup":56,"./Tilelayer":58,"./Tileset":60,"zlibjs":7}],60:[function(_dereq_,module,exports){
var utils = _dereq_('../utils/utils'),
    inherit = _dereq_('../utils/inherit'),
    math = _dereq_('../math/math'),
    Texture = _dereq_('../display/Texture'),
    Vector = _dereq_('../math/Vector'),
    PIXI = _dereq_('pixi.js');

/**
 * This object represents a tileset used by a Tilemap.
 * There can be multiple Tilesets in a map
 *
 * @class Tileset
 * @extends Texture
 * @constructor
 * @param texture {Texture} The texture to use for the tileset
 * @param settings {Object} All the settings for the tileset
 * @param settings.tilewidth {Number} The width of a single tile in the set
 * @param settings.tileheight {Number} The height of a single tile in the set
 * @param [settings.firstgid=1] {Number} The id of the first tile in the set, defaults to 1
 * @param [settings.spacing=0] {Number} The spacing around tiles in the tileset (in pixels)
 * @param [settings.margin=0] {Number} The margin around a tile in the tileset (in pixels)
 * @param [settings.tileoffset] {Object} The offset to apply to a tile rendered from this tileset
 * @param [settings.tileoffset.x=0] {Number} The X offset to apply to the tile
 * @param [settings.tileoffset.y=0] {Number} The Y offset to apply to the tile
 * @param [settings.properties] {Object} User-defined, custom properties that apply to the tileset
 * @param [settings.tileproperties] {Object} User-defined, custom properties that apply to tiles in the tileset.
 *          The keys of this object should the tile id of the properties
 * @param [settings.imagewidth] {Number} An override for the image width
 * @param [settings.imageheight] {Number} An override for the image height
 */
//TODO: Support external tilesets (TSX files) via the "source" attribute
//see: https://github.com/bjorn/tiled/wiki/TMX-Map-Format#tileset
var Tileset = function(texture, settings) {
    //initialize the base Texture class
    if(texture instanceof Array) {
        this.multi = true;
        Texture.call(this, texture[0].baseTexture);
    } else {
        Texture.call(this, texture.baseTexture || texture);
    }

    //Tiled Editor properties

    /**
     * The first tileId in the tileset
     *
     * @property firstgid
     * @type Number
     */
    this.firstgid = settings.firstgid || 1;

    /**
     * The name of the tileset
     *
     * @property name
     * @type String
     */
    this.name = settings.name;

    /**
     * The size of a tile in the tileset
     *
     * @property tileSize
     * @type Vector
     */
    this.tileSize = new Vector(settings.tilewidth, settings.tileheight);

    /**
     * The spacing around a tile in the tileset
     *
     * @property spacing
     * @type Number
     */
    this.spacing = settings.spacing || 0;

    /**
     * The margin around a tile in the tileset
     *
     * @property margin
     * @type Number
     */
    this.margin = settings.margin || 0;

    /**
     * The offset of tile positions when rendered
     *
     * @property tileoffset
     * @type Number
     */
    this.tileoffset = new Vector(
        settings.tileoffset ? settings.tileoffset.x : 0,
        settings.tileoffset ? settings.tileoffset.y : 0
    );

    //TODO: Support for "terraintypes," "image"
    //see: https://github.com/bjorn/tiled/wiki/TMX-Map-Format#tileset

    //Custom/Optional properties

    /**
     * The number of tiles calculated based on size, margin, and spacing
     *
     * @property numTiles
     * @type Vector
     */
    this.numTiles = this.multi ? texture.length : new Vector(
        math.floor((this.baseTexture.source.width - this.margin) / (this.tileSize.x - this.spacing)),
        math.floor((this.baseTexture.source.height - this.margin) / (this.tileSize.y - this.spacing))
    );

    /**
     * The last tileId in the tileset
     *
     * @property lastgid
     * @type Number
     */
    this.lastgid = this.firstgid + (this.multi ? texture.length : ((this.numTiles.x * this.numTiles.y) || 1)) - 1;

    /**
     * The properties of the tileset
     *
     * @property properties
     * @type Object
     */
    this.properties = settings.properties || {};

    /**
     * The properties of the tiles in the tileset (like collision stuff)
     *
     * @property tileproperties
     * @type Object
     */
    this.tileproperties = settings.tileproperties || {};

    /**
     * The size of the tileset
     *
     * @property size
     * @type Vector
     */
    this.size = this.multi ? Vector.ZERO : new Vector(
        settings.imagewidth || this.baseTexture.source.width,
        settings.imageheight || this.baseTexture.source.height
    );

    /**
     * The texture instances for each tile in the set
     *
     * @property textures
     * @type Array
     */
    this.textures = this.multi ? texture : [];

    //massages strings into the values they should be
    //i.e. "true" becomes the value: true
    this.properties = utils.parseTiledProperties(this.properties);

    //massage tile properties
    for(var k in this.tileproperties) {
        this.tileproperties[k] = utils.parseTiledProperties(this.tileproperties[k]);
    }

    //generate tile textures
    if(!this.multi) {
        for(var t = 0, tl = this.lastgid - this.firstgid + 1; t < tl; ++t) {
            //convert the tileId to x,y coords of the tile in the Texture
            var y = math.floor(t / this.numTiles.x),
                x = (t - (y * this.numTiles.x));

            //get location in pixels
            x = (x * this.tileSize.x) + (x * this.spacing) + this.margin;
            y = (y * this.tileSize.y) + (y * this.spacing) + this.margin;

            this.textures.push(
                new Texture(
                    this.baseTexture,
                    new PIXI.Rectangle(x, y, this.tileSize.x, this.tileSize.y)
                )
            );
        }
    }
};

inherit(Tileset, Texture, {
    /**
     * Gets the tile properties for a tile based on it's ID
     *
     * @method getTileProperties
     * @param tileId {Number} The id of the tile to get the properties for
     * @return {Object} The properties of the tile
     */
    getTileProperties: function(tileId) {
        if(!tileId) return null;

        var flags = Tileset.FLAGS,
            flippedX = tileId & flags.FlippedX,
            flippedY = tileId & flags.FlippedY,
            rotatedCW = tileId & flags.RotatedCW;

        tileId = (tileId & ~Tileset.FLAGS.ALL) - this.firstgid;

        //if less than 0, then this id isn't in this tileset
        if(tileId < 0) return null;

        var props = this.tileproperties[tileId] ?
        //get this value
        this.tileproperties[tileId] :
        //set this id to default values and cache
        this.tileproperties[tileId] = {
            collidable: false,
            breakable: false
        };

        props.flippedX = flippedX;
        props.flippedY = flippedY;
        props.rotatedCW = rotatedCW;

        return props;
    },
    /**
     * Gets the tile texture for a tile based on it's ID
     *
     * @method getTileTexture
     * @param tileId {Number} The id of the tile to get the texture for
     * @return {Texture} The texture for the tile
     */
    getTileTexture: function(tileId) {
        if(!tileId) return null;

        //get the internal ID of the tile in this set (0 indexed)
        tileId = (tileId & ~Tileset.FLAGS.ALL) - this.firstgid;

        //if less than 0, then this id isn't in this tileset
        if(tileId < 0) return null;

        return this.textures[tileId];
    },
    /**
     * Returns whether or not this tileset contains the given tile guid
     *
     * @method contains
     * @param tileId {Number} The ID of the tile to check
     * @return {Boolean}
     */
    contains: function(tileId) {
        if(!tileId) return false;

        tileId &= ~Tileset.FLAGS.ALL;

        return (tileId >= this.firstgid && tileId <= this.lastgid);
    }
});

/**
 * Tileset GID flags, these flags are set on a tile's ID to give it a special property
 *
 * @property FLAGS
 * @static
 */
Tileset.FLAGS = {
    FlippedX: 0x80000000,
    FlippedY: 0x40000000,
    RotatedCW: 0x20000000
};

var mask = 0;
for(var f in Tileset.FLAGS) {
    mask |= Tileset.FLAGS[f];
}

Tileset.FLAGS.ALL = mask;

module.exports = Tileset;

},{"../display/Texture":21,"../math/Vector":47,"../math/math":48,"../utils/inherit":69,"../utils/utils":71,"pixi.js":6}],61:[function(_dereq_,module,exports){
// Thanks to PhotonStorm (http://photonstorm.com/) for this loader!
// heavily insprite by (stolen from): https://github.com/photonstorm/phaser/blob/master/src/loader/Cache.js

var inherit = _dereq_('./inherit'),
    C = _dereq_('../constants'),
    Texture = _dereq_('../display/Texture'),
    BaseTexture = _dereq_('../display/BaseTexture'),
    BitmapText = _dereq_('../text/BitmapText'),
    Tilemap = _dereq_('../tilemap/Tilemap'),
    PIXI = _dereq_('pixi.js');

/**
 * A game only has one instance of a Cache and it is used to store all externally loaded assets such
 * as images, sounds and data files as a result of Loader calls. Cache items use string based keys for look-up.
 *
 * @class Cache
 * @extends Object
 * @constructor
 * @param game {Game} The game instance this cache belongs to
 */
var Cache = function(game) {
    /**
     * Local reference to Game.
     *
     * @property game
     * @type Game
     */
    this.game = game;

    /**
     * Canvas key-value container.
     *
     * @property _canvases
     * @type Object
     * @private
     */
    this._canvases = {};

    /**
     * Image key-value container.
     *
     * @property _images
     * @type Object
     * @private
     */
    this._images = {};

    /**
     * Sound key-value container.
     *
     * @property _sounds
     * @type Object
     * @private
     */
    this._sounds = {};

    /**
     * Text key-value container.
     *
     * @property _text
     * @type Object
     * @private
     */
    this._text = {};

    /**
     * Tilemap key-value container.
     *
     * @property _tilemaps
     * @type Object
     * @private
     */
    this._tilemaps = {};

    this.addDefaultImage();
};

inherit(Cache, Object, {
    /**
     * Add a new canvas.
     *
     * @method addCanvas
     * @param obj {Object} The spritesheet object
     * @param obj.key {String} Asset key for this canvas.
     * @param obj.canvas {HTMLCanvasElement} Canvas DOM element.
     * @param obj.context {CanvasRenderingContext2D} Render context of this canvas.
     */
    addCanvas: function(obj) {
        this._canvases[obj.key] = obj;
    },

    /**
     * Add a new sprite sheet.
     *
     * @method addSpriteSheet
     * @param obj {Object} The spritesheet object
     * @param obj.key {String} Asset key for the sprite sheet.
     * @param obj.url {String} URL of this sprite sheet file.
     * @param obj.image {Image} The image of the sprite sheet
     * @param obj.frameWidth {number} Width of the sprite sheet.
     * @param obj.frameHeight {number} Height of the sprite sheet.
     * @param obj.frameMax {number} How many frames stored in the sprite sheet.
     */
    addSpriteSheet: function(obj) {
        var key = obj.key;

        PIXI.BaseTextureCache[key] = new BaseTexture(obj.image);
        PIXI.TextureCache[key] = new Texture(PIXI.BaseTextureCache[key]);
        obj.texture = PIXI.TextureCache[key];

        obj.textures = Texture.fromSpritesheet(obj);

        this._images[key] = obj;
    },

    /**
     * Add a new tilemap.
     *
     * @method addTilemap
     * @param obj {Object} The tilemap file object
     * @param obj.key  {String} Asset key for the tilemap
     * @param obj.url  {String} URL of the tilemap data file
     * @param obj.data {Object} The loaded tilemap data
     * @param obj.format {Number} The format of the tilemap data
     * @param [obj.images] {Array<Image>} Array of images used in the tilesets of this tilemap
     */
    addTilemap: function(obj) {
        //parse out an object representing the XML map
        if(obj.format === C.FILE_FORMAT.XML) {
            obj.xmlData = obj.data;
            obj.data = Tilemap.parseXMLMap(obj.data);
        }

        //create the textures for this map
        var tsets = obj.data.tilesets,
            tset = null,
            name = '',
            k, k2;

        obj.textures = {};

        for(var i = 0, il = tsets.length; i < il; ++i) {
            tset = tsets[i];

            name = tset.name;
            k = obj.key + '_' + name;

            if(tset.image) {

                PIXI.BaseTextureCache[k] = new BaseTexture(obj.images[tset.image]);
                PIXI.TextureCache[k] = new Texture(PIXI.BaseTextureCache[k]);

                obj.textures[name] = PIXI.TextureCache[k];
            } else if(tset.tiles) {
                obj.textures[name] = [];
                for(var t in tset.tiles) {
                    k2 = k + '_' + t;

                    PIXI.BaseTextureCache[k2] = new BaseTexture(obj.images[tset.tiles[t].image]);
                    PIXI.TextureCache[k2] = new Texture(PIXI.BaseTextureCache[k2]);

                    obj.textures[name][t] = PIXI.TextureCache[k2];
                }
            }
        }

        this._tilemaps[obj.key] = obj;
    },

    /**
     * Add a new texture atlas.
     *
     * @method addTextureAtlas
     * @param obj {Object} The texture atlas file object
     * @param obj.key  {String} Asset key for the texture atlas.
     * @param obj.url  {String} URL of this texture atlas file.
     * @param obj.format {Number} The format of the atlas data ATLAS_FORMAT.JSON_ARRAY, ATLAS_FORMAT.JSON_HASH, or ATLAS_FORMAT.STARLING_XML
     * @param obj.data {Object} The texture atlas data exported from TexturePacker
     * @param obj.image {Image} The texture image
     */
    addTextureAtlas: function(obj) {
        var key = obj.key;

        PIXI.BaseTextureCache[key] = new BaseTexture(obj.image);
        PIXI.TextureCache[key] = new Texture(PIXI.BaseTextureCache[key]);
        obj.texture = PIXI.TextureCache[key];

        if(obj.format === C.ATLAS_FORMAT.JSON_ARRAY || obj.format === C.ATLAS_FORMAT.JSON_HASH) {
            obj.textures = Texture.fromJSON(key, obj.data, obj.texture.baseTexture);
        }
        else if (obj.format ===  C.ATLAS_FORMAT.STARLING_XML) {
            obj.textures = Texture.fromXML(key, obj.data, obj.texture.baseTexture);
        }

        this._images[key] = obj;
    },

    /**
     * Add a new Bitmap Font.
     *
     * @method addBitmapFont
     * @param obj {Object} The bitmap font file object
     * @param obj.key  {String} Asset key for the font texture.
     * @param obj.url  {String} URL of this font xml file.
     * @param obj.data {Object} Extra font data.
     * @param obj.format {Number} The format of the bitmap font data
     */
    addBitmapFont: function(obj) {
        var key = obj.key;

        PIXI.BaseTextureCache[key] = new BaseTexture(obj.image);
        PIXI.TextureCache[key] = new Texture(PIXI.BaseTextureCache[key]);
        obj.texture = PIXI.TextureCache[key];

        obj.font = BitmapText.parseXMLFont(key, obj.data, obj.texture);

        this._images[key] = obj;
    },

    /**
     * Add a new image.
     *
     * @method addImage
     * @param obj {Object} The image file object
     * @param obj.key {String} Asset key for the image.
     * @param obj.url {String} URL of this image file.
     * @param obj.image {Image} The image object that was loaded
     */
    addImage: function(obj) {
        var key = obj.key;

        PIXI.BaseTextureCache[key] = new BaseTexture(obj.image);
        PIXI.TextureCache[key] = new Texture(PIXI.BaseTextureCache[key]);
        obj.texture = PIXI.TextureCache[key];

        this._images[key] = obj;
    },

    /**
     * Add a new sound.
     *
     * @method addAudio
     * @param obj {Object} The audio file object
     * @param obj.key {String} Asset key for the audio.
     * @param obj.url {String} URL of this audio file.
     * @param obj.data {ArrayBuffer|Audio} The loaded audio data
     * @param obj.webAudio {Boolean} Is this a webAudio ArrayBuffer for a sound?
     * @param obj.decoded {Boolean} Is the data decoded yet?
     */
    addAudio: function(obj) {
        var key = obj.key;

        if(!obj.webAudio) {
            obj.decoded = true;
        }

        obj.isDecoding = false;

        this._sounds[key] = obj;
    },

    updateSound: function(key, property, value) {
        if(this._sounds[key])
            this._sounds[key][property] = value;
    },

    /**
     * Add a new text data.
     *
     * @method addText
     * @param obj {Object} The text file object
     * @param obj.key {String} Asset key for the text data.
     * @param obj.url {String} URL of this text data file.
     * @param obj.data {object} Extra text data.
     */
    addText: function(obj) {
        this._text[obj.key] = obj;
    },

    /**
     * Adds a default image to be used when a key is wrong / missing.
     * Is mapped to the key __default
     */
    addDefaultImage: function () {
        var key = '__default';

        var base = new BaseTexture();
        base.width = 0;
        base.height = 0;
        base.hasLoaded = true; // avoids a hanging event listener

        PIXI.BaseTextureCache[key] = base;
        PIXI.TextureCache[key] = new Texture(base);

        Texture.__default = PIXI.TextureCache[key];

        this._images[key] = {
            texture: PIXI.TextureCache[key]
        };
    },

    /**
     * Get canvas by key.
     *
     * @method getCanvas
     * @param key {String} Asset key of the canvas you want.
     * @return {HTMLCanvasElement}
     */
    getCanvas: function(key) {
        if(this._canvases[key])
            return this._canvases[key].canvas;
    },

    /**
     * Get image data by key.
     *
     * @method getImage
     * @param key {String} Asset key of the image you want.
     * @return {Image}
     */
    getImage: function(key) {
        if(this._images[key])
            return this._images[key].image;
    },

    /**
     * Get a Texture by key.
     *
     * @method
     * @param key {String} Asset key of the RenderTexture you want.
     * @return {Texture}
     */
    getTexture: function(key) {
        if(this._images[key])
            return this._images[key].texture;
    },

    /**
     * Get a Texture by key.
     *
     * @method
     * @param key {String} Asset key of the RenderTexture you want.
     * @return {Texture}
     */
    getTextures: function(key) {
        if(this._images[key])
            return this._images[key].textures;
    },

    /**
     * Get a Bitmap Font by key.
     *
     * @method
     * @param key {String} Asset key of the Bitmap Font you want.
     * @return {Texture}
     */
    getBitmapFont: function(key) {
        if(this._images[key])
            return this._images[key].font;
    },

    /**
     * Get tilemap data by key.
     *
     * @method getTilemap
     * @param key {String} Asset key of the tilemap you want.
     * @return {object} The tilemap file data. The map data is in the `data` property, the images (for tileset) are in `images`
     */
    getTilemap: function(key) {
        return this._tilemaps[key];
    },

    /**
     * Get sound by key.
     *
     * @method getAudio
     * @param key {String} Asset key of the sound you want.
     * @return {Object}
     */
    getAudio: function(key) {
        return this._sounds[key];
    },

    /**
     * Get sound data by key.
     *
     * @method getAudioData
     * @param key {String} Asset key of the sound you want.
     * @return {ArrayBuffer|Audio}
     */
    getAudioData: function(key) {
        if(this._sounds[key])
            return this._sounds[key].data;
    },

    /**
     * Get text data by key.
     *
     * @method getText
     * @param key {String} Asset key of the text data you want.
     * @return {object} The text data you want.
     */
    getText: function(key) {
        if(this._text[key])
            return this._text[key].data;
    },

    /**
     * Remove a canvas by key.
     *
     * @method removeCanvas
     * @param key {String} key to remove
     */
    removeCanvas: function(key) {
        delete this._canvases[key];
    },

    /**
     * Remove an image by key.
     *
     * @method removeImage
     * @param key {String} key to remove
     */
    removeImage: function(key) {
        delete this._images[key];
    },

    /**
     * Remove a sound by key.
     *
     * @method removeSound
     * @param key {String} key to remove
     */
    removeSound: function(key) {
        delete this._sounds[key];
    },

    /**
     * Remove a text by key.
     *
     * @method removeText
     * @param key {String} key to remove
     */
    removeText: function(key) {
        delete this._text[key];
    },

    /**
     * Destroys this object, removing references so the GC can cleanup
     *
     * @method destroy
     */
    destroy: function() {
        //lose references to let GC cleanup
        this.game = null;
        this._canvases = null;
        this._images = null;
        this._sounds = null;
        this._text = null;
        this._tilemaps = null;
    }
});

module.exports = Cache;

},{"../constants":11,"../display/BaseTexture":15,"../display/Texture":21,"../text/BitmapText":54,"../tilemap/Tilemap":59,"./inherit":69,"pixi.js":6}],62:[function(_dereq_,module,exports){
var inherit = _dereq_('./inherit');

/**
 * High performance clock, based on mrdoob's
 * [Three.js clock](https://github.com/mrdoob/three.js/blob/master/src/core/Clock.js),
 * but with tweaks.
 *
 * @class Clock
 * @extends Object
 * @constructor
 */
var Clock = function() {
    this.startTime = 0;
    this.oldTime = 0;
    this.elapsedTime = 0;

    this.running = false;

    this.timer = window.performance && window.performance.now ? window.performance : Date;
};

inherit(Clock, Object, {
    /**
     * Gets the current time from the underlying timer
     *
     * @method now
     * @return {Number} The current timestamp
     * @example
     *      clock.now();
     */
    now: function() {
        return this.timer.now();
    },
    /**
     * Starts the timer
     *
     * @method start
     * @return {Clock} Returns itself.
     * @chainable
     * @example
     *      clock.start();
     */
    start: function() {
        this.startTime = this.oldTime = this.now();
        this.running = true;

        return this;
    },
    /**
     * Stops the timer
     *
     * @method stop
     * @return {Clock} Returns itself.
     * @chainable
     * @example
     *      clock.stop();
     */
    stop: function() {
        this.getElapsedTime();
        this.running = false;

        return this;
    },
    /**
     * Resets the timer
     *
     * @method reset
     * @return {Clock} Returns itself.
     * @chainable
     * @example
     *      clock.reset();
     */
    reset: function() {
        this.elapsedTime = 0;
        this.startTime = this.oldTime = this.now();

        return this;
    },
    /**
     * Gets the total time that the timer has been running
     *
     * @method getElapsedTime
     * @return {Number} Total ellapsed time in ms
     * @example
     *      clock.getElapsedTime();
     */
    getElapsedTime: function() {
        this.getDelta();

        return this.elapsedTime;
    },
    /**
     * Gets the difference in time since getDelta() was called last
     *
     * @method getDelta
     * @return {Number} Ellapsed time since last call in seconds
     * @example
     *      clock.getDelta();
     */
    getDelta: function() {
        var diff = 0;

        if(this.running) {
            var newTime = this.now();

            diff = 0.001 * (newTime - this.oldTime);
            this.oldTime = newTime;

            this.elapsedTime += diff;
        }

        return diff;
    }
});

module.exports = Clock;

},{"./inherit":69}],63:[function(_dereq_,module,exports){
var inherit = _dereq_('./inherit'),
    math = _dereq_('../math/math');

var Color = function(r, g, b, a) {
    this._color = 0x000000;
    this._rgba = {
        r: 0,
        g: 0,
        b: 0,
        a: 255
    };

    //one arg is either a hex color or an rgba object
    if(arguments.length === 1) {
        //hex color
        if(typeof r === 'number') {
            this.color = r;
        }
        //rgba object
        else {
            this.rgba = r;
        }
    }
    //otherwise real rgba
    else {
        this.setRgba(r, g, b, a);
    }
};

inherit(Color, Object, {
    setRgba: function(r, g, b, a) {
        this._rgba.r = r || 0;
        this._rgba.g = g || 0;
        this._rgba.b = b || 0;
        this._rgba.a = a || 255;

        this._color = (this._rgba.r << 16) + (this._rgba.g << 8) + this._rgba.b;
    }
});

Object.defineProperty(Color.prototype, 'rgba', {
    get: function() {
        return this._rgba;
    },
    set: function(obj) {
        obj = obj || {};
        this._rgba.r = obj.r || 0;
        this._rgba.g = obj.g || 0;
        this._rgba.b = obj.b || 0;
        this._rgba.a = obj.a || 255;

        this._color = (this._rgba.r << 16) + (this._rgba.g << 8) + this._rgba.b;
    }
});

Object.defineProperty(Color.prototype, 'color', {
    get: function() {
        return this._color;
    },
    set: function(clr) {
        this._color = clr;

        this._rgba.r = (clr & 0xff0000) >> 16;
        this._rgba.g = (clr & 0x00ff00) >> 8;
        this._rgba.b = (clr & 0x0000ff);
    }
});

Object.defineProperty(Color.prototype, 'alpha', {
    get: function() {
        return this._rgba.a;
    },
    set: function(alpha) {
        this._rgba.a = math.clamp(alpha, 0, 255);
    }
});

module.exports = Color;

},{"../math/math":48,"./inherit":69}],64:[function(_dereq_,module,exports){
/**
 * Event emitter mixin. This will add emitter properties to an object so that
 * it can emit events, and have others listen for them. Based on
 * [node.js event emitter](https://github.com/joyent/node/blob/master/lib/events.js)
 *
 * @class EventEmitter
 * @constructor
 */
var EventEmitter = function() {
    this._events = this._events || {};

    /**
     * Registers a listener function to be run on an event occurance
     *
     * @method on
     * @param type {String} The event name to listen for
     * @param listener {Function} The function to execute when the event happens
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.addEventListener = this.on = function(type, listener) {
        if(typeof listener !== 'function')
            throw new TypeError('listener must be a function');

        if(!this._events)
            this._events = {};

        // Optimize the case of one listener. Don't need the extra array object.
        if (!this._events[type])
            this._events[type] = listener;
        // If we've already got an array, just append.
        else if (typeof this._events[type] === 'object')
            this._events[type].push(listener);
        // Adding the second element, need to change to array.
        else
            this._events[type] = [this._events[type], listener];

        return this;
    };

    /**
     * Emits an event which will run all registered listeners for the event type
     *
     * @method emit
     * @param type {String} The event name to emit
     * @param data {mixed} Any data you want passed along with the event
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.dispatchEvent = this.emit = function(type) {
        if(!this._events)
            this._events = {};

        var handler = this._events[type],
            len = arguments.length,
            htype = typeof handler,
            args,
            i,
            listeners;

        if(htype === 'undefined')
            return false;

        if(htype === 'function') {
            switch(len) {
                // fast cases
                case 1:
                    handler.call(this);
                    break;
                case 2:
                    handler.call(this, arguments[1]);
                    break;
                case 3:
                    handler.call(this, arguments[1], arguments[2]);
                    break;
                // slower
                default:
                    args = new Array(len - 1);
                    for(i = 1; i < len; ++i)
                        args[i - 1] = arguments[i];

                    handler.apply(this, args);
                    break;
            }
        } else if (htype === 'object') {
            args = new Array(len - 1);
            for(i = 1; i < len; ++i)
                args[i - 1] = arguments[i];

            /*listeners = handler.slice();
            len = listeners.length;
            for (i = 0; i < len; i++)
                listeners[i].apply(this, args);*/
            len = handler.length;
            listeners = handler.slice();
            for(i = 0; i < len; ++i)
                listeners[i].apply(this, args);
        }

        return this;
    };

    /**
     * Removes a listener function for an event type
     *
     * @method off
     * @param type {String} The event name to emit
     * @param listener {Function} The function to remove
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.removeEventListener = this.off = function(type, listener) {
        var list, position, length, i;

        if(typeof listener !== 'function')
            throw new TypeError('listener must be a function');

        if(!this._events[type])
            return this;

        list = this._events[type];
        length = list.length;
        position = -1;

        if(list === listener || (typeof list.listener === 'function' && list.listener === listener)) {
            delete this._events[type];
        } else if(typeof list === 'object') {
            for(i = length; i-- > 0;) {
                if(list[i] === listener || (list[i].listener && list[i].listener === listener)) {
                    position = i;
                    break;
                }
            }

            if(position < 0)
                return this;

            if(list.length === 1) {
                list.length = 0;
                delete this._events[type];
            } else {
                list.splice(position, 1);
            }
        }

        return this;
    };

    /**
     * Registers a one-time callback for an event
     *
     * @method once
     * @param type {String} The event name to listen for
     * @param listener {Function} the callback to call when the event occurs
     * @return {mixed} Returns itself.
     * @chainable
     */
    this.once = function(type, listener) {
        if(typeof listener !== 'function')
            throw new TypeError('listener must be a function');

        var fired = false;

        function g() {
            this.off(type, g);

            if(!fired) {
                fired = true;
                listener.apply(this, arguments);
            }
        }

        g.listener = listener;
        this.on(type, g);

        return this;
    };
};

module.exports = EventEmitter;

},{}],65:[function(_dereq_,module,exports){
var inherit = _dereq_('./inherit'),
    Sprite = _dereq_('../display/Sprite'),
    Tilemap = _dereq_('../tilemap/Tilemap'),
    Rectangle = _dereq_('../geom/Rectangle'),
    BitmapText = _dereq_('../text/BitmapText');

/**
 * The object factory makes it simple to create and add objects to a parent. One is added
 * to a State's world and camera by default, but they can be used for any parent but they
 * can only belong to a single state.
 *
 * @class ObjectFactory
 * @extends Object
 * @constructor
 * @param state {State} The game state this factory belongs to
 * @param parent {Container} The container to act as the parent for created objects
 */
var ObjectFactory = function(state, parent) {
    this.state = state;
    this.game = state.game;
    this.parent = parent;
};

inherit(ObjectFactory, Object, {
    /**
     * Adds a generic object to the world or camera
     *
     * @method obj
     * @param object {mixed} Any game object you want to add to the parent
     * @return {mixed} Returns the added object
     */
    obj: function(obj) {
        return this.parent.addChild(obj);
    },
    /**
     * Creates a new sprite and adds it to the game world
     *
     * @method sprite
     * @param texture {String|Texture} The texture for the sprite, or the key for one in the cache
     * @param [frame=null] {String|Number} A specific frame of a sprite sheet to use, either the index or string key
     *      depending on the type of the sheet when loaded.
     * @param [physics=true] {Boolean} Should this sprite be added to the physics simulation?
     * @return {Sprite} The sprite added
     */
    sprite: function(tx, frame, physics) {
        var spr,
            game = this.game;

        if(typeof tx === 'string') {
            if(frame || frame === 0)
                tx = game.cache.getTextures(tx)[frame];
            else
                tx = game.cache.getTexture(tx);
        }

        if(!tx) {
            tx = game.cache.getTexture('__default');
        }

        spr = new Sprite(tx);

        //if undefined, then default to true
        if(physics || physics === undefined) {
            spr.enablePhysics(this.state.physics);
            //this.state.physics.addSprite(spr);
        }

        return this.parent.addChild(spr);
    },
    /**
     * Creates a new AudioPlayer to play the sound passed in
     *
     * @method audio
     * @param key {String} The unique cache key for the preloaded audio
     * @param [settings] {Object} All the settings for the audio player (see AudioManager.add for all settings)
     * @return {AudioPlayer} The player added
     */
    audio: function(key, settings) {
        return this.state.audio.add(key, settings);
    },
    /**
     * Creates a new tilemap to add to the world
     *
     * @method tilemap
     * @param key {String} The unique cache key for the preloaded tilemap data
     * @param [constrain=true] {Boolean} Should the camera be constrained to this tilemap's size?
     * @return {Tilemap} The tilemap added
     */
    tilemap: function(key, constrain) {
        var obj = this.game.cache.getTilemap(key) || {},
            tilemap = new Tilemap(this.state, obj.data, obj.textures);

        if(constrain) {
            this.state.camera.constrain(new Rectangle(0, 0, tilemap.realSize.x, tilemap.realSize.y));
        }

        //force render of tilemap
        tilemap.render(
            -this.state.world.position.x,
            -this.state.world.position.x,
            this.game.width,
            this.game.height
        );

        tilemap._cachekey = key;

        return this.parent.addChild(tilemap);
    },
    /**
     * Creates a new instance of BitmapText
     *
     * @method bitmaptext
     * @param text {String} The text for the BitmapText to display
     * @param font {String} The key for the bitmap font loaded into the cache
     * @param interactive {Boolean} Can the item be interacted with by mouse (clicked, dragged, etc)
     * @return {BitmapText} The bitmap text object added
     */
    bitmaptext: function(text, font, style) {
        if(typeof font === 'string')
            font = this.game.cache.getBitmapFont(font);

        return this.parent.addChild(new BitmapText(text, font, style));
    }
});

module.exports = ObjectFactory;

},{"../display/Sprite":19,"../geom/Rectangle":37,"../text/BitmapText":54,"../tilemap/Tilemap":59,"./inherit":69}],66:[function(_dereq_,module,exports){
var inherit = _dereq_('./inherit');

/**
 * Holds a pool of different Objects to help reduce the number times
 * an object is created and destroyed.
 *
 * @class ObjectPool
 * @extends Object
 * @constructor
 * @param type {mixed} The object type that this pool will hold (like Sprite, or Tile)
 * @param parent {mixed} The parent that the objects will be added to. Passing this in will
 *      make the pool add any newly created objects as children to this object.
 */
var ObjectPool = function(type, parent) {
    this.type = type;
    this.pool = [];
    this.parent = parent;
};

inherit(ObjectPool, Object, {
    /**
     * Creates a new instance of the pool's object type, or if available
     * pulls one that is already created out of the pool
     *
     * @method create
     * @return {mixed} The instance of the object pulled from the pool
     */
    create: function() {
        var o = this.pool.pop();

        if(!o) {
            o = this._construct(this.type, arguments);
            if(this.parent)
                this.parent.addChild(o);
        }

        o.__allocated = true;

        return o;
    },
    /**
     * Frees an object back into the pool to be recycled
     *
     * @method free
     */
    free: function(o) {
        //don't free something twice
        if(o.__allocated) {
            o.__allocated = false;
            this.pool.push(o);
        }
    },
    //have to do this hack around to be able to use
    //apply and new together
    _construct: function(ctor, args) {
        function F() {
            return ctor.apply(this, args);
        }
        F.prototype = ctor.prototype;
        return new F();
    }
});

module.exports = ObjectPool;

},{"./inherit":69}],67:[function(_dereq_,module,exports){
var inherit = _dereq_('./inherit');

var Queue = function() {
    this.data = [];
    this.index = 0;
};

inherit(Queue, Object, {
    enqueue: function(item) {
        this.data.push(item);
    },

    dequeue: function() {
        if(this.data.length === 0)
            return;

        var item = this.data[this.index];

        if(++this.index * 2 >= this.data.length) {
            this.data = this.data.slice(this.index);
            this.index = 0;
        }

        return item;
    },

    peek: function() {
        return (this.data.length > 0 ? this.data[this.index] : undefined);
    }
});

Object.defineProperty(Queue.prototype, 'length', {
    get: function() {
        return (this.data.length - this.index);
    }
});

Object.defineProperty(Queue.prototype, 'empty', {
    get: function() {
        return (this.data.length === 0);
    }
});

module.exports = Queue;

},{"./inherit":69}],68:[function(_dereq_,module,exports){
var inherit = _dereq_('./inherit'),
    Sprite = _dereq_('../display/Sprite');

/**
 * Holds a pool of different Sprites that can be created, makes it very
 * easy to quickly create different registered entities
 *
 * @class SpritePool
 * @extends Object
 * @constructor
 * @param game {Game} The game instance this sprite pool belongs to
 */
var SpritePool = function(game) {
    this.types = {};

    /**
     * The game instance this pool belongs to
     *
     * @property game
     * @type Game
     */
    this.game = game;

    this.add('_default', Sprite);
};

inherit(SpritePool, Object, {
    /**
     * Adds an Sprite Type to the pool
     *
     * @method add
     * @param name {String} The user-defined name of the Sprite Type to add
     * @param obj {Sprite} The Sprite or decendant type to add to the pool
     * @return {Sprite} Returns the passed sprite
     */
    add: function(name, obj) {
        return this.types[name] = obj;
    },
    /**
     * Checks if the Sprite Type exists in the pool
     *
     * @method has
     * @param name {String} The user-defined name of the Sprite Type to check if is in the pool
     * @return {Boolean}
     */
    has: function(name) {
        return !!this.types[name];
    },
    /**
     * Creates a new sprite from the pool
     *
     * @method create
     * @param name {String} The user-defined name of the Sprite to check if is in the pool
     * @param texture {Texture} The texture for the sprite
     * @param props {Object} Extra object that will be passed along (for custom sprite options)
     * @return {Sprite} Returns a new instance of the object from the pool
     */
    create: function(name, texture, props) {
        if(!name || !this.types[name])
            name = '_default';

        return new this.types[name](texture, props);
    },
    /**
     * Frees an object back into the pool to be recycled, currently doesn't do any recycling unfortunately
     *
     * @method free
     */
    free: function() {
        return;
    }
});

module.exports = SpritePool;

},{"../display/Sprite":19,"./inherit":69}],69:[function(_dereq_,module,exports){
/**
 * Inherits the prototype of a parent object.
 *
 * @method inherit
 * @param child {Function} The Child to inherit the prototype
 * @param parent {Function} The Parent to inherit from
 * @param proto {Object} The prototype to apply to the child
 */
var inherit = function(child, parent, proto) {
    proto = proto || {};

    //get the property descriptors from the child proto and the passed proto
    var desc = {};
    [child.prototype, proto].forEach(function (s) {
        Object.getOwnPropertyNames(s).forEach(function (k) {
            desc[k] = Object.getOwnPropertyDescriptor(s, k);
        });
    });

    //set the constructor descriptor
    desc.constructor = {
        value: child,
        enumerable: false,
        writable: true,
        configurable: true
    };

    //create the prototype
    child.prototype = Object.create(parent.prototype, desc);
};

module.exports = inherit;
},{}],70:[function(_dereq_,module,exports){
/**
 * Holds the results of the feature detection run on the browser, to make it simple to
 * see which features the library can use.
 *
 * @class support
 * @extends Object
 * @static
 */
var support = {
    /**
     * The current user agent string
     *
     * @property ua
     * @type String
     */
    ua: window.navigator ? window.navigator.userAgent.toLowerCase() : 'nodejs',

    /**
     * Whether or not canvas is supported
     *
     * @property canvas
     * @type Boolean
     */
    canvas: !!(function () { try { return window.CanvasRenderingContext2D && document.createElement('canvas').getContext('2d'); } catch(e) { return false; } })(),

    /**
     * Whether or not webgl is supported
     *
     * @property webgl
     * @type Boolean
     */
    webgl: !!(function () { try { var c = document.createElement('canvas'); return window.WebGLRenderingContext && (c.getContext('webgl') || c.getContext('experimental-webgl')); } catch(e) { return false; } })(),

    /**
     * Whether or not the crypto API is supported
     *
     * @property crypto
     * @type Boolean
     */
    crypto: !!window.crypto && !!window.crypto.getRandomValues,

    /**
     * Whether or not web workers are supported
     *
     * @property workers
     * @type Boolean
     */
    workers: !!window.Worker,

    /**
     * Whether or not Blob URLs are supported
     *
     * @property blobs
     * @type Boolean
     */
    blobUrls: !!window.Blob && !!window.URL && !!window.URL.createObjectURL,

    /**
     * Whether or not typed arrays are supported
     *
     * @property typedArrays
     * @type Boolean
     */
    typedArrays: !!window.ArrayBuffer,

    /**
     * Whether or not the filesystem API is supported
     *
     * @property fileapi
     * @type Boolean
     */
    fileapi: !!window.File && !!window.FileReader && !!window.FileList && !!window.Blob,

    /**
     * Whether or not the Web Audio API is supported
     *
     * @property webAudio
     * @type Boolean
     */
    webAudio: !!window.AudioContext || !!window.webkitAudioContext || !!window.mozAudioContext,

    /**
     * Whether html Audio is supported in this browser
     *
     * @property htmlAudio
     * @type Boolean
     */
    htmlAudio: !!document.createElement('audio').canPlayType && !!window.Audio,

    /**
     * Whether or not local storage is supported
     *
     * @property localStorage
     * @type Boolean
     */
    localStorage: !!window.localStorage,

    /**
     * Whether or not touch is supported
     *
     * @property touch
     * @type Boolean
     */
    touch: !!(('createTouch' in document) || ('ontouchstart' in window) || (navigator.isCocoonJS)),

    /**
     * Whether or not the gamepad API is supported
     *
     * @property gamepad
     * @type Boolean
     */
    gamepad: !!navigator.webkitGetGamepads || !!navigator.webkitGamepads || (navigator.userAgent.indexOf('Firefox/') !== -1)
};

/**
 * Describes which audio codecs a browser supports
 *
 * @property codec
 * @type Object
 */
if(support.htmlAudio) {
    var audioTest = new Audio();

    support.codec = {
        mp3: !!audioTest.canPlayType('audio/mpeg;').replace(/^no$/,''),
        opus: !!audioTest.canPlayType('audio/ogg; codecs="opus"').replace(/^no$/,''),
        ogg: !!audioTest.canPlayType('audio/ogg; codecs="vorbis"').replace(/^no$/,''),
        wav: !!audioTest.canPlayType('audio/wav; codecs="1"').replace(/^no$/,''),
        m4a: !!(audioTest.canPlayType('audio/x-m4a;') || audioTest.canPlayType('audio/aac;')).replace(/^no$/,''),
        webm: !!audioTest.canPlayType('audio/webm; codecs="vorbis"').replace(/^no$/,'')
    };
}

module.exports = support;

},{}],71:[function(_dereq_,module,exports){
var Vector = _dereq_('../math/Vector'),
    Circle = _dereq_('../geom/Circle'),
    Rectangle = _dereq_('../geom/Rectangle'),
    Polygon = _dereq_('../geom/Polygon');

/**
 * The grapefruit utility object, used for misc functions used throughout the code base
 *
 * @class utils
 * @extends Object
 * @static
 */
var utils = {
    _arrayDelim: /[|,]/,
    /**
     * An empty function that performs no action
     *
     * @method noop
     */
    noop: function() {},
    /**
     * Gets the absolute url from a relative one
     *
     * @method getAbsoluteUrl
     * @param url {String} The relative url to translate into absolute
     * @return {String} The absolute url (fully qualified)
     */
    getAbsoluteUrl: function(url) {
        var a = document.createElement('a');
        a.href = url;
        return a.href;
    },
    /**
     * Performs an ajax request, and manages the callbacks passed in
     *
     * @method ajax
     * @param settings {Object} The settings of the ajax request, similar to jQuery's ajax function
     * @return {XMLHttpRequest|ActiveXObject} An XHR object
     */
    ajax: function(sets) {
        //base settings
        sets = sets || {};
        sets.method = sets.method || 'GET';
        sets.dataType = sets.dataType || 'text';

        if(!sets.url)
            throw new TypeError('Undefined URL passed to ajax');

        //callbacks
        sets.progress = sets.progress || utils.noop;
        sets.load = sets.load || utils.noop;
        sets.error = sets.error || utils.noop;
        sets.abort = sets.abort || utils.noop;
        sets.complete = sets.complete || utils.noop;

        var xhr = utils.createAjaxRequest(),
            protocol = utils.getAbsoluteUrl(sets.url).split('/')[0];

        xhr.onreadystatechange = function() {
            if(xhr.readyState === 4) {
                var res = xhr.response || xhr.responseText,
                    err = null;

                //The 'file:' protocol doesn't give response codes
                if(protocol !== 'file:' && xhr.status !== 200)
                    err = 'Non-200 status code returned: ' + xhr.status;

                if(!err && typeof res === 'string') {
                    if(sets.dataType === 'json') {
                        try {
                            res = JSON.parse(res);
                        } catch(e) {
                            err = e;
                        }
                    } else if(sets.dataType === 'xml') {
                        try {
                            res = utils.parseXML(res);
                        } catch(e) {
                            err = e;
                        }
                    }
                }

                if(err) {
                    if(sets.error) sets.error.call(xhr, err);
                } else {
                    if(sets.load) sets.load.call(xhr, res);
                }
            }
        };

        xhr.open(sets.method, sets.url, true);

        //chrome doesn't support json responseType, some browsers choke on XML type
        if(sets.dataType === 'arraybuffer')
            xhr.responseType = 'arraybuffer';
        else
            xhr.responseType = 'text';

        xhr.send();

        return xhr;
    },
    /**
     * Wraps XMLHttpRequest in a cross-browser way.
     *
     * @method AjaxRequest
     * @return {XMLHttpRequest|ActiveXObject}
     */
    //from pixi.js
    createAjaxRequest: function() {
        //activeX versions to check for in IE
        var activexmodes = ['Msxml2.XMLHTTP', 'Microsoft.XMLHTTP'];

        //Test for support for ActiveXObject in IE first (as XMLHttpRequest in IE7 is broken)
        if(window.ActiveXObject) {
            for(var i=0; i<activexmodes.length; i++) {
                try {
                    return new window.ActiveXObject(activexmodes[i]);
                }
                catch(e) {
                    //suppress error
                }
            }
        }
        // if Mozilla, Safari etc
        else if(window.XMLHttpRequest) {
            return new XMLHttpRequest();
        }
        else {
            return false;
        }
    },
    /**
     * This will take values and override the passed obj's properties with those values.
     * The difference from a normal object extend is that this will try to massage the passed
     * value into the same type as the object's property. Also if the key for the value is not
     * in the original object, it is not copied.
     *
     * @method setValues
     * @param obj {Object} The object to extend the values into
     * @param values {Object} The values to put into the object
     * @return {Object} returns the updated object
     * @example
     *      var obj = { vec: new Vector(), arr: [] },
     *          vals = { vec: '2|5', arr: '5|10|11' };
     *      utils.setValues(obj, vals);
     *      //now obj is:
     *      // { vec: Vector(2, 5), arr: [5, 10, 11] }
     *
     */
    //similar to https://github.com/mrdoob/three.js/blob/master/src/materials/Material.js#L42
    setValues: function(obj, values) {
        if(!values) return;

        for(var key in values) {
            var newVal = values[key];

            if(newVal === undefined) {
                //console.warn('Object parameter '' + key + '' is undefined.');
                continue;
            }
            if(key in obj) {
                var curVal = obj[key];

                //massage strings into numbers
                if(typeof curVal === 'number' && typeof newVal === 'string') {
                    var n;
                    if(newVal.indexOf('0x') === 0) n = parseInt(newVal, 16);
                    else n = parseInt(newVal, 10);

                    if(!isNaN(n))
                        obj[key] = n;
                    /*else
                        console.warn('Object parameter '' + key + '' evaluated to NaN, using default. Value passed: ' + newVal);*/

                }
                //massage vectors
                else if(curVal instanceof Vector && newVal instanceof Array) {
                    curVal.set(parseFloat(newVal[0], 10) || 0, parseFloat(newVal[1], 10) || parseFloat(newVal[0], 10) || 0);
                } else if(curVal instanceof Vector && typeof newVal === 'string') {
                    var a = newVal.split(utils._arrayDelim, 2);
                    curVal.set(parseFloat(a[0], 10) || 0, parseFloat(a[1], 10) || parseFloat(a[0], 10) || 0);
                } else if(curVal instanceof Vector && typeof newVal === 'number') {
                    curVal.set(newVal, newVal);
                }
                //massage arrays
                else if(curVal instanceof Array && typeof newVal === 'string') {
                    obj[key] = newVal.split(utils._arrayDelim);
                    for(var i = 0, il = obj[key].length; i < il; ++i) {
                        var val = obj[key][i];
                        if(!isNaN(val)) obj[key][i] = parseFloat(val, 10);
                    }
                } else {
                    obj[key] = newVal;
                }
            }
        }

        return obj;
    },
    /**
     * From jQuery.extend, extends one object into another
     * taken straight from jQuery 2.0.3
     *
     * @method extend
     */
    extend: function() {
        var src, copyIsArray, copy, name, options, clone, target = arguments[0] || {},
            i = 1,
            length = arguments.length,
            deep = false;

        // Handle a deep copy situation
        if (typeof target === 'boolean') {
            deep = target;
            target = arguments[1] || {};
            // skip the boolean and the target
            i = 2;
        }

        // Handle case when target is a string or something (possible in deep copy)
        if (typeof target !== 'object' && typeof target !== 'function') {
            target = {};
        }

        // extend jQuery itself if only one argument is passed
        //if (length === i) {
        //    target = this;
        //    --i;
        //}

        for (; i < length; i++) {
            // Only deal with non-null/undefined values
            options = arguments[i];
            if (options !== null && options !== undefined) {
                // Extend the base object
                for (name in options) {
                    src = target[name];
                    copy = options[name];

                    // Prevent never-ending loop
                    if (target === copy) {
                        continue;
                    }

                    // Recurse if we're merging plain objects or arrays
                    if (deep && copy && (utils.isPlainObject(copy) || (copyIsArray = Array.isArray(copy)))) {
                        if (copyIsArray) {
                            copyIsArray = false;
                            clone = src && Array.isArray(src) ? src : [];

                        } else {
                            clone = src && utils.isPlainObject(src) ? src : {};
                        }

                        // Never move original objects, clone them
                        target[name] = utils.extend(deep, clone, copy);

                        // Don't bring in undefined values
                    } else if (copy !== undefined) {
                        target[name] = copy;
                    }
                }
            }
        }

        // Return the modified object
        return target;
    },
    /**
     * From jQuery.isPlainObject, checks if an object is a plain object
     * taken straight from jQuery 2.0.3
     *
     * @method isPlainObject
     * @param obj {mixed} The object to test
     * @return {Boolean}
     */
    isPlainObject: function(obj) {
        // Not plain objects:
        // - Any object or value whose internal [[Class]] property is not "[object Object]"
        // - DOM nodes
        // - window
        if (typeof obj !== 'object' || obj.nodeType || obj === obj.window) {
            return false;
        }

        // Support: Firefox <20
        // The try/catch suppresses exceptions thrown when attempting to access
        // the "constructor" property of certain host objects, ie. |window.location|
        // https://bugzilla.mozilla.org/show_bug.cgi?id=814622
        try {
            if (obj.constructor && !Object.hasOwnProperty.call(obj.constructor.prototype, 'isPrototypeOf')) {
                return false;
            }
        } catch(e) {
            return false;
        }

        // If the function hasn't returned already, we're confident that
        // |obj| is a plain object, created by {} or constructed with new Object
        return true;
    },
    /**
     * Get the DOM offset values of any given element
     *
     * @method getOffset
     * @param element {HTMLElement} The targeted element that we want to retrieve the offset
     * @return {Vector} The offset of the element
     */
    getOffset: function(element) {
        var box = element.getBoundingClientRect(),
            clientTop = element.clientTop || document.body.clientTop || 0,
            clientLeft = element.clientLeft || document.body.clientLeft || 0,
            scrollTop = window.pageYOffset || element.scrollTop || document.body.scrollTop,
            scrollLeft = window.pageXOffset || element.scrollLeft || document.body.scrollLeft;

        return new Vector(
            box.left + scrollLeft - clientLeft,
            box.top + scrollTop - clientTop
        );
    },
    /**
     * Parses an array of numbers that represent a hitArea into the actual shape.
     *
     * For example: `[1, 1, 15]` is a Circle (`[x, y, radius]`); `[1, 1, 15, 15]` is a Rectangle
     * (`[x, y, width, height]`); and any length >= 5 is a polygon in the form `[x1, y1, x2, y2, ..., xN, yN]`.
     *
     * @method parseHitArea
     * @param value {Array<Number>} The array to parse
     * @return {Circle|Rectangle|Polygon} The parsed out shape
     */
    parseHitArea: function(hv) {
        var ha;

        //odd number of values
        if(hv.length % 2 !== 0 && hv.length !== 3) {
            throw new RangeError('Strange number of values for hitArea! Should be a flat array of values, like: [x,y,r] for a circle, [x,y,w,h] for a rectangle, or [x,y,x,y,...] for other polygons.');
        }

        //a circle x,y,r
        if(hv.length === 3) {
            ha = new Circle(hv[0], hv[1], hv[2]);
        }
        //a rectangle x,y,w,h
        else if(hv.length === 4) {
            ha = new Rectangle(hv[0], hv[1], hv[2], hv[3]);
        }
        //generic polygon
        else {
            ha = new Polygon(0, 0, hv);
        }

        return ha;
    },
    /**
     * Parses an object of string properties into potential javascript types. First it attempts to
     * convert to a number, if that fails it checks for the string 'true' or 'false' and changes it
     * to the actual Boolean value, then it attempts to parse a string as JSON.
     *
     * @method parseTiledProperties
     * @param value {Array<Number>} The array to parse
     * @return {Circle|Rectangle|Polygon} The parsed out shape
     */
    parseTiledProperties: function(obj) {
        if(!obj || obj.__tiledparsed)
            return obj;

        for(var k in obj) {
            var v = obj[k],
                n = parseFloat(v, 10);

            //try to massage numbers
            if(n === 0 || n)
                obj[k] = n;
            //true values
            else if(v === 'true')
                obj[k] = true;
            //false values
            else if(v === 'false')
                obj[k] = false;
            //anything else is either a string or json, try json
            else {
                try{
                    v = JSON.parse(v);
                    obj[k] = v;
                } catch(e) {}
            }
        }

        //after parsing, check some other things
        if(obj.hitArea)
            obj.hitArea = utils.parseHitArea(obj.hitArea);

        if(obj.body === 'static' || obj.sensor) {
            obj.mass = Infinity;
            obj.inertia = Infinity;
        }

        obj.__tiledparsed = true;

        return obj;
    },
    _logger: window.console || {},
    /**
     * Safe way to log to console, if console.log doesn't exist nothing happens.
     *
     * @method log
     */
    log: function() {
        if(utils._logger.log)
            utils._logger.log.apply(utils._logger, arguments);
    },
    /**
     * Safe way to warn to console, if console.warn doesn't exist nothing happens.
     *
     * @method warn
     */
    warn: function() {
        if(utils._logger.warn)
            utils._logger.warn.apply(utils._logger, arguments);
    },
    /**
     * Safe way to error to console, if console.error doesn't exist nothing happens.
     *
     * @method error
     */
    error: function() {
        if(utils._logger.error)
            utils._logger.error.apply(utils._logger, arguments);
    }
};

/**
 * Parses an XML string into a Document object. Will use window.DOMParser
 * if available, falling back to Microsoft.XMLDOM ActiveXObject in IE.
 *
 * Eventually, it would be nice to include a node.js alternative as well
 * for running in that environment.
 *
 * @method parseXML
 * @param xmlStr {String} The xml string to parse
 * @return {Document} An XML Document
 */
//XML Parser in window
if(typeof window.DOMParser !== 'undefined') {
    utils.parseXML = function(xmlStr) {
        return (new window.DOMParser()).parseFromString(xmlStr, 'text/xml');
    };
}
//IE specific XML parser
else if(typeof window.ActiveXObject !== 'undefined' && new window.ActiveXObject('Microsoft.XMLDOM')) {
    utils.parseXML = function(xmlStr) {
        var xmlDoc = new window.ActiveXObject('Microsoft.XMLDOM');
        xmlDoc.async = 'false';
        xmlDoc.loadXML(xmlStr);
        return xmlDoc;
    };
}
//node.js environment
/*else if(__isNode) {
    utils.parseXML = function(xmlStr) {
        var DOMParser = require('xmldom').DOMParser;
        return (new DOMParser()).parseFromString(xmlStr, "text/xml");
    };
}*/
// no parser available
else {
    utils.warn('XML parser not available, trying to parse any XML will result in an error.');
    utils.parseXML = function() {
        throw new Error('Trying to parse XML, but not XML parser is available in this environment');
    };
}

module.exports = utils;

},{"../geom/Circle":34,"../geom/Polygon":36,"../geom/Rectangle":37,"../math/Vector":47}]},{},[14])
(14)
});