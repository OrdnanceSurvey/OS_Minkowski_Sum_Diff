/*
 * This is an experimental implementation of Minkowski sum and difference based on JTS geometric functionality.
 * 
 * Current implementation supports Minkowski sums between a "source" geometry (polygon/linestring/multipolygon/multilinestring/GeometryCollection) and a "reference" geometry (polygon/linestring), 
 * and Minkowski difference between a "source" geometry (polygon/multipolygon/GeometryCollection) and a "reference" geometry (polygon/linestring) 
 * 
 * Polygons may be concave. The "source" polygon may contain holes. 
 * 
 * Any holes in "reference" polygon are ignored (in most cases it doesn't make practical sense anyway).
 * 
 * You will need JTS 1.15 to use this implementation. It will NOT work with JTS1.14 or earlier due to JTS package name issues
 *
 * version 0.3
 * 
 * date: 2019-01-19
 * 
 * author: Sheng Zhou (Sheng.Zhou@os.uk)
 * 
 * Copyright (C) 2019 Ordnance Survey
 *
 * Licensed under the Open Government Licence v3.0 (the "License");
 * 
 * you may not use this file except in compliance with the License.
 * 
 * You may obtain a copy of the License at
 *
 *     http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */
//
package uk.osgb.algorithm.minkowski_sum;

import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Collection;
import java.util.Vector;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;

public class Minkowski_Sum {
	//
	public static GeometryFactory gf = new GeometryFactory();
	//
	/*************************************************************
	 * 
	 *  Minkowski Sum of a Polygon/MultiPolygon/LineString/MultiLineString and a Polygon/LineString
	 * 
	 *************************************************************/
	/** Minkowski sum of two general geometries, dispatched by manually testing types 
	 * @param src source geometry which might be polygon/multipolygon/linestring/multilinestring
	 * @param ref reference geometry, which might be a polygon or a linestring
	 * @param doReflection If reflection over origin is performed first on reference polygon (e.g. for collision detection purpose). 
	 * @return the sum or an empty polygon (for un-supported types)
	 */
	public static Geometry compMinkSum(Geometry src, Geometry ref, boolean doReflection, boolean isRefConvex) {
		if(src==null || ref == null) {
			return gf.createPolygon(); 
		}
		if(ref instanceof Polygon) {
			if(src instanceof Polygon) {
				return compMinkSumPlgPlg((Polygon)src, (Polygon)ref, doReflection, isRefConvex);
			}else if(src instanceof LineString) {
				return compMinkSumLSPlg((LineString)src, (Polygon)ref, doReflection, isRefConvex);
			}else if(src instanceof MultiPolygon) {
				return compMinkSumMultiPlgPlg((MultiPolygon)src, (Polygon)ref, doReflection, isRefConvex);
			}else if(src instanceof MultiLineString) {
				return compMinkSumMultiLSPlg((MultiLineString)src, (Polygon)ref, doReflection, isRefConvex);
			}else if(src instanceof GeometryCollection) {
				return compMinkSumGeometryCollection((GeometryCollection) src, ref, doReflection, isRefConvex);
			}else if(src instanceof Point) {
				return compMinkSumPoint((Point) src, ref, doReflection, isRefConvex);
			}
		}else if(ref instanceof LineString || ref instanceof LinearRing) {
			if(src instanceof Polygon) {
				return compMinkSumPlgLS((Polygon)src, (LineString)ref, doReflection);
			}else if(src instanceof LineString) {
				return compMinkSumLSLS((LineString)src, (LineString)ref, doReflection);
			}else if(src instanceof MultiPolygon) {
				return compMinkSumMultiPlgLS((MultiPolygon)src, (LineString)ref, doReflection);
			}else if(src instanceof MultiLineString) {
				return compMinkSumMultiLSLS((MultiLineString)src, (LineString)ref, doReflection);
			}else if(src instanceof GeometryCollection) {
				return compMinkSumGeometryCollection((GeometryCollection) src, ref, doReflection, isRefConvex);
			}else if(src instanceof Point) {
				return compMinkSumPoint((Point) src, ref, doReflection, isRefConvex);
			}
		}
		return src.getFactory().createPolygon();
	}
    //
	//trying out simulated multi-dispatch
	//
	/** Minkowski sum of two general geometries, using simulated multi-dispatch
	 * @param src source geometry which might be polygon/multipolygon/linestring/multilinestring
	 * @param ref reference geometry, which might be a polygon or a linestring
	 * @param doReflection If reflection over origin is performed first on reference polygon (e.g. for collision detection purpose). 
	 * @param isRefConvex if reference geometry is known to be convex, a small bit of exter performance improvement may be achieved
	 * @return minkowski sum of src and ref (or an empty geometry)
	 */
	public static Geometry minkSum(Geometry src, Geometry ref, Boolean doReflection, Boolean isRefConvex) {
		if(src==null || ref == null) {
			return gf.createPolygon(); 
		}
		try {
			Class[] clzs = new Class[] {
					src.getClass(), ref.getClass(), Boolean.class, Boolean.class
			};
			Object[] objs = new Object[] {
					src, ref, doReflection, isRefConvex
			};
			Class clz = Class.forName("uk.osgb.algorithm.minkowski_sum.Minkowski_Sum"); 
			Method m = clz.getDeclaredMethod("minkSum", clzs);

			return (Geometry) m.invoke(clz, objs);

		} catch (IllegalAccessException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return src.getFactory().createPolygon(); // return an empty geometry
	}
	/**  Minkowski sum of two general geometries, computed without any presumption
	 * @param src source geometry which might be polygon/multipolygon/linestring/multilinestring
	 * @param ref reference geometry, which might be a polygon or a linestring
	 * @return minkowski sum of src and ref (or an empty geometry)
	 */
	public static Geometry minkSum(Geometry src, Geometry ref) {
		return minkSum(src, ref, false, false);
	}
	// wrappers for testing
	private static Geometry minkSum(Point src, Polygon ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumPoint(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(Point src, LineString ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumPoint(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(Point src, LinearRing ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumPoint(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(Polygon src, Polygon ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumPlgPlg(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(LineString src, Polygon ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumLSPlg(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(MultiPolygon src, Polygon ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumMultiPlgPlg(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(MultiLineString src, Polygon ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumMultiLSPlg(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(Polygon src, LineString ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumPlgLS(src, ref, doReflection);
	}
	private static Geometry minkSum(Polygon src, LinearRing ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumPlgLS(src, ref, doReflection);
	}
	private static Geometry minkSum(LineString src, LineString ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumLSLS(src, ref, doReflection);
	}
	private static Geometry minkSum(LineString src, LinearRing ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumLSLS(src, ref, doReflection);
	}
	private static Geometry minkSum(MultiPolygon src, LineString ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumMultiPlgLS(src, ref, doReflection);
	}
	private static Geometry minkSum(MultiPolygon src, LinearRing ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumMultiPlgLS(src, ref, doReflection);
	}
	private static Geometry minkSum(MultiLineString src, LineString ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumMultiLSLS(src, ref, doReflection);
	}
	private static Geometry minkSum(MultiLineString src, LinearRing ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumMultiLSLS(src, ref, doReflection);
	}
	private static Geometry minkSum(GeometryCollection src, Polygon ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumGeometryCollection(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(GeometryCollection src, LineString ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumGeometryCollection(src, ref, doReflection, isRefConvex);
	}
	private static Geometry minkSum(GeometryCollection src, LinearRing ref, Boolean doReflection, Boolean isRefConvex) {
		return compMinkSumGeometryCollection(src, ref, doReflection, isRefConvex);
	}
	//
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumPoint(Point src, Geometry ref, boolean doReflection, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric on coordinate origin, the reflected geometry is the same as original
			if(ref instanceof Polygon) {
				ref = reflectionOrgPlg((Polygon) ref);
			}else if(ref instanceof LineString || ref instanceof LinearRing) {
				ref = reflectionOrgLS((LineString) ref);
			}
		}
		AffineTransformation af = new AffineTransformation();
		af.translate(src.getX(), src.getY());
		return af.transform(ref);
	}
	
	public static Geometry compMinkSumMultiPoint(MultiPoint src, Geometry ref, boolean doReflection, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric on coordinate origin, the reflected geometry is the same as original
			if(ref instanceof Polygon) {
				ref = reflectionOrgPlg((Polygon) ref);
			}else if(ref instanceof LineString || ref instanceof LinearRing) {
				ref = reflectionOrgLS((LineString) ref);
			}
		}
		Geometry sum = src.getFactory().createPolygon();
		int numPt = src.getNumGeometries();
		for(int i = 0; i < numPt; ++i) {
			Point pt = (Point) src.getGeometryN(i);
			AffineTransformation af = new AffineTransformation();
			af.translate(pt.getX(), pt.getY());
			Geometry geom = af.transform(ref);
			sum = sum.union(geom);
		}
		return sum;
	}
	/** Minkowski sum of two polygons. Any holes in reference polygon are ignored
	 * @param src Source polygon, may contain holes
	 * @param ref Reference polygon, holes are ignored
	 * @param doReflection If reflection is performed first on reference polygon (e.g. for collision detection purpose).
	 * For repeated computation, reflection may be performed first before calling this method.
	 * @param isRefConvex whether the reference polygon is convex (if unknown, please use false)
	 * @return
	 */
	public static Geometry compMinkSumPlgPlg(Polygon src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric on coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkSumPlg(src, refCoords, isRefConvex);
	}
	//
	/** Minkowski sum of a multipolygon and a polygon
	 * @param src source multipolygon, may contain holes
	 * @param ref reference polygon, holes are ignored
	 * @param doReflection if reflection is performed first on reference polygon (for repeated computation, reflection may be performed first before calling this method)
	 * @param isRefConvex whether the reference polygon is convex (if unknown, please use false)
	 * @return
	 */
	public static Geometry compMinkSumMultiPlgPlg(MultiPolygon src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgPlg(ref);
		}
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry sum1 = compMinkSumPlgPlg(plg, ref, false, isRefConvex);
			try {
				rlt = rlt.union(sum1);
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}
	/** Minkowski sum between a polygon and a reference linestring
	 * @param src Soource polygon, may contain holes
	 * @param ref reference linestring
	 * @param doReflection if reflection is performed first on reference linestring (e.g. for collision detection purpose).
	 * For repeated computation, reflection may be performed first before calling this method.
	 * @return
	 */
	public static Geometry compMinkSumPlgLS(Polygon src, LineString ref, boolean doReflection) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgLS(ref);
		}
		//
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkSumPlg(src, refCoords, false); // always treat LS ref as non-convex
	}
	//
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	public static Geometry compMinkSumMultiPlgLS(MultiPolygon src, LineString ref, boolean doReflection) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric over coordinate origin, the  geometry is the same as original
			ref = reflectionOrgLS(ref);
		}
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry sum1 = compMinkSumPlgLS(plg, ref, false);
			try {
				rlt = rlt.union(sum1);
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumLSPlg(LineString src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkSumLS(src, refCoords, isRefConvex);
	}
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	public static Geometry compMinkSumLSLS(LineString src, LineString ref, boolean doReflection) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkSumLS(src, refCoords, false);

	}//
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	public static Geometry compMinkSumMultiLSLS(MultiLineString src, LineString ref, boolean doReflection) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			LineString ls = (LineString) src.getGeometryN(i);
			Geometry sum1 = compMinkSumLS(ls, refCoords, false);
			try {
				rlt = rlt.union(sum1);
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}
	
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumMultiLSPlg(MultiLineString src, Polygon ref, boolean doReflection, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doReflection) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			LineString ls = (LineString) src.getGeometryN(i);
			Geometry sum1 = compMinkSumLS(ls, refCoords, isRefConvex);
			try {
				rlt = rlt.union(sum1);
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}
	/** Minkowskie sum of an arbitrary GeometryCollection and a polygon
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumGeometryCollection(GeometryCollection src, Geometry ref, boolean doReflection, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Geometry part =  src.getGeometryN(i);
			Geometry sum1 = minkSum(part, ref, doReflection, isRefConvex);
			try {
				rlt = rlt.union(sum1);
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
		return rlt;
	}
	//
	
	//
	/** Minkowski sum between a linestring and a reference array of coordinates, assuming src is not null or emtpy (handled elsewehere)
	 * @param src source linestring
	 * @param refCoords coordinate sequence, which forms either an open LineString or a closed LinearRing
	 * @param isRefConvex if refCoords is a closed LinearRing, whether it is convex or not
	 * @return
	 */
	private static Geometry compMinkSumLS(LineString src, Coordinate[] refCoords, boolean isRefConvex) {
		Coordinate[] coords = src.getCoordinates();
		// expanding the shell of original geometry (on two directions)
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring
			return coordArrayVectorAddition(coords, refCoords, isRefConvex, src.getFactory());
		}else {
			return coordArrayVectorAddition(coords, refCoords, src.getFactory());
		}
	}
	/**
	 * @param src
	 * @param refCoords
	 * @param isRefConvex 
	 * @return
	 */
	private static Geometry compMinkSumPlg(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if(src==null||src.isEmpty()) {
			return src;
		}
		boolean isRefRing = false;
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring
			isRefRing = true;
		}
		//
		LineString extRing = src.getExteriorRing();
		Coordinate[] extCoords = extRing.getCoordinates();
		Geometry extSum = expansionShell(gf.createPolygon(extCoords), refCoords, isRefConvex);
		int numHoles = src.getNumInteriorRing();		
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				LineString holeRing = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeRing.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords);
				Geometry intDiff = shrinkHole(holePlg, refCoords, isRefConvex);
				extSum = extSum.difference(intDiff);
			}
		}
		return extSum;
	}

	/**Compute the "vector addition" of refCoords and the exterior ring of src and return the exterior ring of the result (including holes outside original src geometry)
	 * @param src Polygon
	 * @param refCoords 
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry expansionShell(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if(src==null||src.isEmpty()) {
			return src;
		}
		boolean isRefRing = false;
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring
			isRefRing = true;
		}
		//
		Geometry sum = null;
		//
		LineString extRing = src.getExteriorRing();
		Coordinate[] extCoords = extRing.getCoordinates();
		// expanding the shell of original geometry (on two directions)
		if(isRefRing) {
			sum = coordArrayVectorAddition(extCoords, refCoords, isRefConvex, src.getFactory());
		}else {
			sum = coordArrayVectorAddition(extCoords, refCoords, src.getFactory());
		}
		// if holes are outside original geometry, they are generated by boundary expansion and should be kept as holes in final results
		Vector<LinearRing> ringVec = null;

		LinearRing[] shellHoles = getHoles((Polygon)sum);
		if(shellHoles!=null) {
			Polygon geomNoHole = gf.createPolygon(extCoords);// polygon formed by the exterior ring of src geometry
			for(int i = 0; i < shellHoles.length; ++i) {
				//Point ep = shellHoles[i].getEndPoint(); // a point on the hole boundary
				Point ep = gf.createPolygon(shellHoles[i]).getInteriorPoint();
				if(!geomNoHole.contains(ep)){// holes outside original geometry, generated by expansion, should be part of the result as a hole
				//if(!geomNoHole.intersects(ep)){// holes outside original geometry, generated by expansion, should be part of the result as a hole
					if(ringVec == null) {
						ringVec = new Vector<LinearRing>();
					}
					ringVec.add(shellHoles[i]);
				}
			}
		}
		// assemble the result
		LinearRing shell = getShellRing((Polygon)sum);
		if(ringVec!=null) {
			// holes of sumHoles
			LinearRing[] holes = new LinearRing[ringVec.size()];
			ringVec.toArray(holes);
			return gf.createPolygon(shell, holes);
		}else {// no holes as result of expansion from original holes
			return gf.createPolygon(shell.getCoordinates());
		}		
	} 
	//
	/******************************************************************
	 * 
	 *  Minkowski Difference between Polygon/MultiPolygon and Polyon/LineString
	 * 
	 * 
	 *****************************************************************/
    //
	//trying out simulated multi-dispatch. Pity Java does not support it internally:( still has to support all subclasses manually
	// if src is not empty and ref is emtpy, the difference is src; if src is empty and ref is not empty, the difference is an empty set
	// Miknowski 
	/**
	 * @param src
	 * @param ref
	 * @param isRefSymmetric if the reference geometry is symmetric over origin, it is not necessary to do the reflection
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry minkDiff(Geometry src, Geometry ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		if(src==null) {
			return gf.createPolygon();
		}else if(ref==null) {
			return src;
		}
		try {
			Class[] clzs = new Class[] {
					src.getClass(), ref.getClass(), Boolean.class, Boolean.class
			};
			Object[] objs = new Object[] {
					src, ref, isRefSymmetric, isRefConvex
			};
			Class clz = Class.forName("uk.osgb.algorithm.minkowski_sum.Minkowski_Sum"); 
			Method m = clz.getDeclaredMethod("minkDiff", clzs);

			return (Geometry) m.invoke(clz, objs);
		} catch (IllegalAccessException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		} catch (IllegalArgumentException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		}
		//return src.getFactory().createPolygon(null, null);
		return src.getFactory().createPolygon();
	}
	//
	public static Geometry minkDiff(Geometry src, Geometry ref) {
		return minkDiff(src, ref, false, false);
	}
	//
	private static Geometry minkDiff(Point src, Polygon ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return src.getFactory().createPolygon();	
	}
	private static Geometry minkDiff(MultiPoint src, Polygon ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return src.getFactory().createPolygon();	
	}
	private static Geometry minkDiff(Polygon src, Polygon ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffPlgPlg(src, ref, isRefSymmetric, isRefConvex);	
	}
	private static Geometry minkDiff(MultiPolygon src, Polygon ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffMultiPlgPlg(src, ref, isRefSymmetric, isRefConvex);	
	}
	private static Geometry minkDiff(Polygon src, LineString ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffPlgLS(src, ref, isRefSymmetric);
	}
	private static Geometry minkDiff(Polygon src, LinearRing ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffPlgLS(src, ref, isRefSymmetric);
	}
	private static Geometry minkDiff(MultiPolygon src, LineString ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffMultiPlgLS(src, ref, isRefSymmetric);
	}
	private static Geometry minkDiff(MultiPolygon src, LinearRing ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffMultiPlgLS(src, ref, isRefSymmetric);
	}
	private static Geometry minkDiff(GeometryCollection src, Polygon ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffGeometryCollection(src, ref, isRefSymmetric);
	}
	private static Geometry minkDiff(GeometryCollection src, LineString ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffGeometryCollection(src, ref, isRefSymmetric);
	}
	private static Geometry minkDiff(GeometryCollection src, LinearRing ref, Boolean isRefSymmetric, Boolean isRefConvex) {
		return compMinkDiffGeometryCollection(src, ref, isRefSymmetric);
	}
	/** Minkowski difference of two geometries. 
	 * @param src Source geometry, may be polygon or multipolygon
	 * @param ref Reference geometry, may be polygon or linestring/linearring
	 * @param refSymmetric true if ref is symmetric over origin
	 * @param isRefConvex
	 * @return null (not implemented yet)
	 */
	public static Geometry compMinkDiff(Geometry src, Geometry ref, boolean isRefSymmetric, boolean isRefConvex) {
		if(src==null) {
			return gf.createPolygon();
		}else if(ref==null) {
			return src;
		}
		if(ref instanceof Polygon) {
			if(src instanceof Polygon) {
				return compMinkDiffPlgPlg((Polygon)src, (Polygon)ref, isRefSymmetric, isRefConvex);
			}else if(src instanceof MultiPolygon) {
				return compMinkDiffMultiPlgPlg((MultiPolygon)src, (Polygon)ref, isRefSymmetric, isRefConvex);
			}else if(src instanceof GeometryCollection){
				return compMinkDiffGeometryCollection((GeometryCollection)src, ref, isRefSymmetric);
			}else {
				return src.getFactory().createPolygon();
			}
		}else if(ref instanceof LineString || ref instanceof LinearRing) {
			if(src instanceof Polygon) {
				return compMinkDiffPlgLS((Polygon)src, (LineString)ref, isRefSymmetric);
			}else if(src instanceof MultiPolygon) {
				return compMinkDiffMultiPlgLS((MultiPolygon)src, (LineString)ref, isRefSymmetric);
			}else if(src instanceof GeometryCollection){
				return compMinkDiffGeometryCollection((GeometryCollection)src, ref, isRefSymmetric);
			}else {
				return src.getFactory().createPolygon(); // unsupported source type
			}
		}else {
			return src.getFactory().createPolygon();
		}
	}
	//
	/** Minkowski difference of polygon (may contain holes) and reference geometry, using difference between minkDiff(exterior ring) and minkSum(holes)
	 * 
	 * @param src
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkDiffPlg(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if(src.isEmpty()) {
			return src;
		}
		boolean isRefRing = false;
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring so isRefConvex is considered
			isRefRing = true;
		}
	// bufferring exterior ring of source geometry
		LineString extRing = src.getExteriorRing();
		Geometry extDiff = shrinkHole(gf.createPolygon(extRing.getCoordinates()), refCoords, isRefConvex);
		int numHoles = src.getNumInteriorRing();		
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				LineString holeSrc = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeSrc.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords); // hole polygon
				Geometry intSum = expansionShell(holePlg, refCoords, isRefConvex);
				try {
					extDiff = extDiff.difference(intSum);
				}catch(Exception e) {
					e.printStackTrace();
				}
			}
		}
		return extDiff;
	}

	/**  
	 * @param src polygon (holes, if any, are ignored)
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry shrinkHole(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		if(src.isEmpty()) {
			return src;
		}
		GeometryFactory gf = src.getFactory();
		boolean isRefRing = false;
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring so isRefConvex is considered
			isRefRing = true;
		}
		Geometry sum = null;
		LineString extRing = src.getExteriorRing();
		Coordinate[] extCoords = extRing.getCoordinates();
		if(isRefRing) {
			sum = coordArrayVectorAddition(extCoords, refCoords, isRefConvex, src.getFactory());
		}else {
			sum = coordArrayVectorAddition(extCoords, refCoords, src.getFactory());
		}
		// holes in buffering result
		LinearRing[] shellHoles = getHoles((Polygon)sum);
		if(shellHoles==null) {// no holes, shrinked completely
			return src.getFactory().createPolygon();
		}
		//
	// need to check if the holes are inside the shell of the original geometry
	// holes could be generated from outwards expansion as well
		Vector<LinearRing> ringVec = new Vector<LinearRing>(shellHoles.length); 
		Polygon geomNoHole = gf.createPolygon(extCoords); // original geom 
		for(int i = 0; i < shellHoles.length; ++i) {
			//Point ep = shellHoles[i].getEndPoint();
			Point ep = gf.createPolygon(shellHoles[i]).getInteriorPoint();
			if(geomNoHole.contains(ep)){ // the hole is inside the original source geometry, not a hole generated by outwards expansion
			//if(geomNoHole.intersects(ep)){ // the hole is inside the original source geometry, not a hole generated by outwards expansion
				ringVec.add(shellHoles[i]);
			}
		}
		//
		if(ringVec.size() == 0) {
			return gf.createPolygon();
		}else {
			shellHoles = new LinearRing[ringVec.size()];
			ringVec.toArray(shellHoles);
		}
		// shell of result 
		Geometry rltShell = null;
		if(shellHoles.length > 1) {// multipolygon
			Polygon[] plgs = new Polygon[shellHoles.length];
			for(int i = 0; i < plgs.length; ++i) {
				plgs[i] = gf.createPolygon(shellHoles[i]);
			}
			rltShell = gf.createMultiPolygon(plgs);
		}else {// polygon
			rltShell = gf.createPolygon(shellHoles[0]);
		}
		// holes
		if(rltShell!=null) {
			return rltShell;
		}else {
			return gf.createPolygon();
		}
	}
	/** Minkowski difference of two polygons (src - ref). Holes in ref are ignored. 
	 * @param src Source polygon 
	 * @param ref reference polgyon (without hole)
	 * @param doReflection whether ref should be reflected (for repeated calls, ref may be reflected before calling this method). If ref is symmetric with respect to the coordinate origin, there is no need to perform reflection either.
	 * @PARAM isRefConvex If you know the ref polygon is convex, set this to true will improve performance and robustness.
	 * @return
	 */
	public static Geometry compMinkDiffPlgPlg(Polygon src, Polygon ref, boolean isRefSymmetric, boolean isRefConvex) {
		if(!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkDiffPlg(src, refCoords, isRefConvex);
	}
	
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	public static Geometry compMinkDiffPlgLS(Polygon src, LineString ref, boolean isRefSymmetric) {
		if(!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkDiffPlg(src, refCoords, false);
	}
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkDiffMultiPlgPlg(MultiPolygon src, Polygon ref, boolean isRefSymmetric, boolean isRefConvex) {
		if(!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgPlg(ref);
		}
		Geometry rlt = null;
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry diff1 = compMinkDiffPlgPlg(plg, ref, false, isRefConvex);
			if(rlt == null) {
				rlt = diff1;
			}else {
				try {
					rlt = rlt.union(diff1);
				}catch(Exception e) {
					e.printStackTrace();
				}
			}
		}
		return rlt;
	}
	//
	/**
	 * @param src
	 * @param ref
	 * @param doReflection
	 * @return
	 */
	public static Geometry compMinkDiffMultiPlgLS(MultiPolygon src, LineString ref, boolean isRefSymmetric) {
		if(!isRefSymmetric) {// for geometry symmetric over coordinate origin, the reflected geometry is the same as original
			ref = reflectionOrgLS(ref);
		}
		Geometry rlt = null;
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Polygon plg = (Polygon) src.getGeometryN(i);
			Geometry diff1 = compMinkDiffPlgLS(plg, ref, false);
			if(rlt == null) {
				rlt = diff1;
			}else {
				try {
					rlt = rlt.union(diff1);
				}catch(Exception e) {
					e.printStackTrace();
				}
			}
		}
		return rlt;
	}
	//
	public static Geometry compMinkDiffGeometryCollection(GeometryCollection src, Geometry ref, boolean isRefSymmetric) {
		Geometry rlt = null;
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Geometry geom = src.getGeometryN(i);
			Geometry diff1 = minkDiff(geom, ref, isRefSymmetric, false);
			if(rlt == null) {
				rlt = diff1;
			}else {
				try{
					rlt = rlt.union(diff1);
				}catch(Exception e) {
					e.printStackTrace();
				}
			}
		}
		return rlt;
	}
	//
	/***************************************************************
	 * 
	 *  auxiliary methods
	 * 
	 * 
	 ***************************************************************/
	//
	/** Get all holes in plg as LinearRing and put them into holeCol.
	 * @param plg
	 * @param holeCol
	 * @return
	 */
	public static int getHoles(Polygon plg, Collection<LinearRing> holeCol) {
		int numHoles = plg.getNumInteriorRing();
		if(numHoles > 0) {
			for(int i = 0; i < numHoles; ++i) {
				LinearRing ring = plg.getFactory().createLinearRing(plg.getInteriorRingN(i).getCoordinates());
				holeCol.add(ring);
			}
		}
		return numHoles;
	}
	//
	/**
	 * @param plg
	 * @return
	 */
	public static LinearRing getShellRing(Polygon plg) {
		return plg.getFactory().createLinearRing(plg.getExteriorRing().getCoordinates());
	}
	//
	/**
	 * @param plg
	 * @return
	 */
	public static LinearRing[] getHoles(Polygon plg) {
		int numHoles = plg.getNumInteriorRing();
		if(numHoles>0) {
			LinearRing[] holes = new LinearRing[numHoles];
			for(int i = 0; i < numHoles; ++i) {
				holes[i] = plg.getFactory().createLinearRing(plg.getInteriorRingN(i).getCoordinates());
			}
			return holes;
		}else {
			return null;
		}
	}
	//
	/**
	 * @param mulPlg
	 * @return
	 */
	public static LinearRing[] getHoles(MultiPolygon mulPlg) {
		int numPlg = mulPlg.getNumGeometries();
		if(numPlg>0) {
			Vector<LinearRing> holes = new Vector<LinearRing>();
			for(int i = 0; i < numPlg; ++i) {
				Polygon plg = (Polygon) mulPlg.getGeometryN(i);
				getHoles(plg, holes);
			}
			if(holes.size() > 0) {
				LinearRing[] rlt = new LinearRing[holes.size()];
				holes.toArray(rlt);
				return rlt;
			}else {
				return null;
			}
		}
		return null;
	}
	//
	/**
	 * @param mulPlg
	 * @return
	 */
	public static LinearRing[] getShellRings(MultiPolygon mulPlg) {
		int numPlg = mulPlg.getNumGeometries();
		if(numPlg>0) {
			LinearRing[] shells = new LinearRing[numPlg];
			for(int i = 0; i < numPlg; ++i) {
				Polygon plg = (Polygon) mulPlg.getGeometryN(i);
				shells[i] = mulPlg.getFactory().createLinearRing(plg.getExteriorRing().getCoordinates());
			}
			return shells;
		}else {
			return null;
		}
	}
	/** treat geom as a set of vectors P, generate -P and return the geometry constructed from -P 
	 * @param geom the geometry (Polygon, LineString or LinearRing only at present)
	 * @return reflected geometry (null if type is not supported yet)
	 */
	public static Geometry reflectionOrgGeom(Geometry geom) {
		if(geom instanceof Polygon) {
			return reflectionOrgPlg((Polygon)geom);
		}else if(geom instanceof LineString || geom instanceof LinearRing) {
			return reflectionOrgLS((LineString)geom);
		}
		return null;
	}
	//
	/** Symmetry of JTS polygon with respect to coordinate origin
	 * @param plg
	 * @return
	 */
	public static Polygon reflectionOrgPlg(Polygon plg){
		Coordinate[] extNeg = reflectionOrgCoordArray(plg.getExteriorRing().getCoordinates());
		GeometryFactory gf = plg.getFactory();
		int numHoles = plg.getNumInteriorRing();
		if(numHoles > 0) {
			LinearRing[] holes = new LinearRing[numHoles];
			for(int i = 0; i < numHoles; ++i) {
				Coordinate[] holeNeg = reflectionOrgCoordArray(plg.getInteriorRingN(i).getCoordinates());
				holes[i] = gf.createLinearRing(holeNeg);
			}
			return gf.createPolygon(gf.createLinearRing(extNeg), holes);
		}else {
			return gf.createPolygon(extNeg);
		}
	}
	public static LineString reflectionOrgLS(LineString ls){
		Coordinate[] coordsNeg = reflectionOrgCoordArray(ls.getCoordinates());
		return ls.getFactory().createLineString(coordsNeg);
	}
	/** symmetry of a Coordinate array with respect to coordinate origin
	 * @param coords
	 * @return
	 */
	public static Coordinate[] reflectionOrgCoordArray(Coordinate[] coords) {
		if(coords!=null && coords.length > 0) {
			Coordinate[] coordsNeg = new Coordinate[coords.length];
			for(int i = 0; i < coords.length; ++i) {
				Coordinate coord = coords[i];
				coordsNeg[i] = new Coordinate(-coord.x, -coord.y);
			}
			return coordsNeg;
		}else {
			return null;
		}
	}
	/** Given a line segment segSp-segEp, and coordinate array of a polygon refCoords, return the point set union of:
	 * 		vector sum of polygon refCoords and segSp;
	 * 		vector sum of polygon refCoords and segEp;
	 * 		all points covered by the polygon when it moves from segSp to segEp
	 * 
	 * @param refCoords Coordinate array of the reference polygon (without holes)
	 * @param segSp starting point of the line segment 
	 * @param segEp ending point of the line segment
	 * @param isRefConvex if the reference polygon refCoords is convex. If reference polygon is convex, the computation is simplier and robust
	 * @return
	 */
	public static Geometry segPlgAddition(Coordinate[] refCoords, Coordinate segSp, Coordinate segEp, boolean isRefConvex, GeometryFactory gf) {
		if(isRefConvex) {
			int numCoord = refCoords.length;
			Coordinate[] coords = new Coordinate[numCoord*2];
			for(int i = 0; i < refCoords.length;++i) {
				Coordinate ref = refCoords[i];
				coords[i] = new Coordinate(ref.x+segSp.x, ref.y+segSp.y);
				coords[i+numCoord] = new Coordinate(ref.x+segEp.x, ref.y + segEp.y);
			}
			ConvexHull ch = new ConvexHull(coords, gf);
			return ch.getConvexHull();
		} else {
			int numCoord = refCoords.length;
			Coordinate[] coords = new Coordinate[numCoord*2];
			Coordinate[] coords1 = new Coordinate[numCoord];
			Coordinate[] coords2 = new Coordinate[numCoord];
			for(int i = 0; i < refCoords.length;++i) {
				Coordinate ref = refCoords[i];
				coords[i] = new Coordinate(ref.x+segSp.x, ref.y+segSp.y);
				coords[i+numCoord] = new Coordinate(ref.x+segEp.x, ref.y + segEp.y);
				coords1[i] = (Coordinate) coords[i].clone();
				coords2[i] = (Coordinate) coords[i + numCoord].clone();
			}
			ConvexHull chAll = new ConvexHull(coords, gf);
			ConvexHull ch1 = new ConvexHull(coords1, gf);
			ConvexHull ch2 = new ConvexHull(coords2, gf);
			Geometry core = chAll.getConvexHull().difference(ch1.getConvexHull()).difference(ch2.getConvexHull());
			ConvexHull coreHull = new ConvexHull(core);
			Geometry coreNew = coreHull.getConvexHull();
			Polygon plg1 = gf.createPolygon(coords1);
			Polygon plg2 = gf.createPolygon(coords2);
			return coreNew.union(plg1).union(plg2);
		}
	}
	//
	//
	/** vector addition of a segment and a linestring (open, although it will handle linearing as well)
	 * @param refCoords Coordinate array of the reference linestring (whether it is a ring is checked outside this method)
	 * @param segSp
	 * @param segEp
	 * @return the sum (may be an emtpy polygon)
	 */
	public static Geometry segLSAddition(Coordinate segSp, Coordinate segEp, Coordinate[] refCoords, GeometryFactory gf) {
		int numPts = refCoords.length;
		Geometry sum = gf.createPolygon();
		for(int i = 0; i < numPts-1; ++i) {
			Coordinate refSp = refCoords[i];
			Coordinate refEp = refCoords[i+1];
			Geometry sum1 = segVectorAddition(segSp, segEp, refSp, refEp, gf); // 
			if(sum1!=null) {
				try {
					sum = sum.union(sum1);
				}catch(Exception e) {
					e.printStackTrace();
				}
			}
		}
		return sum;
	}
	/** vector addition of a segment on reference linestring and a segment on the source geometry. This is the basic operation upon which all other methods are built
	 *  Convexhull is used to construct result to improve robustness (avoiding self-intersection, hopefully).
	 * @param segSp
	 * @param segEp
	 * @param refSp
	 * @param refEp
	 * @return a polygon (empty polygon if the sum is 0D or 1D)
	 */
	private static Geometry segVectorAddition(Coordinate segSp, Coordinate segEp, Coordinate refSp, Coordinate refEp, GeometryFactory gf){
		Coordinate[] coords = new Coordinate[4];
		coords[0] = new Coordinate(refSp.x + segSp.x, refSp.y + segSp.y);
		coords[1] = new Coordinate(refEp.x + segSp.x, refEp.y + segSp.y);
		coords[2] = new Coordinate(refEp.x + segEp.x, refEp.y + segEp.y);
		coords[3] = new Coordinate(refSp.x + segEp.x, refSp.y + segEp.y);
		ConvexHull ch = new ConvexHull(coords, gf);
		Geometry hull = ch.getConvexHull();
		if(hull.getDimension() > 1) {
			return hull; // a polygon
		}else {
			return gf.createPolygon(); // return an empty polygon
		}
	}
	/** computing the vector sum of a linestring/polygon (without holes) and a polygon in form of coordinate array.
	 * @param geomCoords Coordinate array of the source linestring/polygon (without holes)
	 * @param refPlgCoords the reference polygon
	 * @param isRefConvex whether the reference polygon is convex
	 * @return
	 */
	private static Geometry coordArrayVectorAddition(Coordinate[] geomCoords, Coordinate[] refPlgCoords, boolean isRefConvex, GeometryFactory gf) {
		Geometry sumAll = gf.createPolygon();
		for(int j = 0; j < geomCoords.length-1; ++j) {
			Coordinate segSp = geomCoords[j];
			Coordinate segEp = geomCoords[j+1];
			if(!segSp.equals2D(segEp)) {
				Geometry sum1 = segPlgAddition(refPlgCoords, segSp, segEp, isRefConvex, gf);
				if(sum1!=null) {
					try {
						sumAll = sumAll.union(sum1);
					}catch(Exception e) {
						e.printStackTrace();
					}
				}
			}
		}
		return sumAll;
	}
	
	/**computing the vector sum of a linestring or polygon (without holes) and a linestring (maybe closed) in form of coordinate array.
	 * @param geomCoords Coordinate array of source geometry (linestring, or polygon without holes)
	 * @param refCoords Coordinate array of reference geometry, an open linestring or a closed lienarring
	 * @return
	 */
	public static Geometry coordArrayVectorAddition(Coordinate[] geomCoords, Coordinate[] refCoords, GeometryFactory gf) {
		int numPts = refCoords.length;
		if(refCoords[0].equals2D(refCoords[numPts-1])) { // linear ring
			return coordArrayVectorAddition(geomCoords, refCoords, false, gf);
		}
		Geometry sumAll = gf.createPolygon();
		for(int j = 0; j < geomCoords.length-1; ++j) {
			Coordinate segSp = geomCoords[j];
			Coordinate segEp = geomCoords[j+1];
			if(!segSp.equals2D(segEp)) {
				Geometry sum1 = segLSAddition(segSp, segEp, refCoords, gf);
				if(sum1!=null) {
					try {
						sumAll = sumAll.union(sum1);
					}catch(Exception e) {
						e.printStackTrace();
					}
				}
			}
		}
		return sumAll;	
	}
	/**
	 * @param original
	 * @return
	 */
	public static Coordinate[] deepCloneCoordArray(Coordinate[] original) {
		if(original != null && original.length > 0) {
			Coordinate[] dup = new Coordinate[original.length];
			for(int i = 0; i < original.length;++i) {
				dup[i] = (Coordinate) original[i].clone();
			}
			return dup;
		}
		return null;
	}
	/**
	 * @return the geometry factory used to build geometry object in Minkowski sum/difference computation. 
	 */
	public static GeometryFactory getGeometryFactory() {
		return gf;
	}
	//
	public static void setGeometryFactory(GeometryFactory geomFactory) {
		gf = geomFactory;
	}
	
	public static void main(String[] args) {
		String wktStr1 = "POLYGON((5 5, 25 5, 25 15, 5 15, 5 5), (10 8, 20 8, 20 12, 10 12, 10 8))";
		String wktStr2 = "POLYGON((-1 -1, 1 -1, 2 0, 1 1, -1 1, -1 -1))";
		String wktStr3 = "POLYGON((5 5, 25 5, 25 15, 5 15, 5 5), (8 8, 13 8, 13 12, 8 12, 8 8), (17 8, 22 8, 22 12, 17 12, 17 8))";
		String wktStr4 = "POLYGON((0 -1, 1 0, 0 1, 0 -1))";
		WKTReader rd = new WKTReader();
		try {
			Geometry geom1 = rd.read(wktStr1);
			Geometry ref00 = rd.read(wktStr2);
			System.out.println("Test 01");
			System.out.println("- A -\n"+geom1.toText());
			System.out.println("- B -\n"+ref00.toText());
			System.out.println("--- A + B ---");
			Geometry sum = minkSum(geom1, ref00, false, false);
			System.out.println(sum.toText());
			System.out.println("--- A + (-B) ---");
			sum = minkSum(geom1, ref00, true, false);
			System.out.println(sum.toText());
			System.out.println("--- A - B ---");
			Geometry diff = minkDiff(geom1, ref00, false, false);
			System.out.println(diff.toText());
			//
			geom1 = rd.read(wktStr3);
			ref00 = rd.read(wktStr4);
			System.out.println("\nTest 02");
			System.out.println("- A -\n"+geom1.toText());
			System.out.println("- B -\n"+ref00.toText());
			System.out.println("--- A + B ---");
			sum = minkSum(geom1, ref00, false, false);
			System.out.println(sum.toText());
			System.out.println("--- A + (-B) ---");
			sum = minkSum(geom1, ref00, true, false);
			System.out.println(sum.toText());
			System.out.println("--- A - B ---");
			diff = minkDiff(geom1, ref00, false, false);
			System.out.println(diff.toText());

		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
