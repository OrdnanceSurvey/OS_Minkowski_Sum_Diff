/* this is an experimental implementation of Minkowski sum and difference based on JTS geometric functionality.
 * 
 * Current implementation supports Minkowski sums between polygon/linestring/multipolygon/multilinestring/GeometryCollection and polygon/linestring, and
 * minkowski difference between polygon/multipolygon/GeometryCollection and polygon/linestring 
 * 
 * It doesn't require polygons to be convex. The "source" polygon may contain holes. 
 * 
 * Any holes in "reference" polygon are ignored (in most cases it doesn't make much practical sense anyway).
 * 
 * You will need JTS 1.15 to use this library.
 *
 * version 0.2
 * 
 * date: 2018-08-08
 * 
 * author: Sheng Zhou (Sheng.Zhou@os.uk)
 * 
 */
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
	/** Minkowski sum of two general geometries. 
	 * @param src source geometry which might be polygon/multipolygon/linestring/multilinestring
	 * @param ref reference geometry, which might be a polygon or a linestring
	 * @param doNeg If negation is performed first on reference polygon (e.g. for collision detection purpose).
	 * @return the sum or null (for un-supported types)
	 */
	public static Geometry compMinkSum(Geometry src, Geometry ref, boolean doNeg, boolean isRefConvex) {
		if(src==null || ref == null) {
			return gf.createPolygon(); 
		}
		if(ref instanceof Polygon) {
			if(src instanceof Polygon) {
				return compMinkSumPlgPlg((Polygon)src, (Polygon)ref, doNeg, isRefConvex);
			}else if(src instanceof LineString) {
				return compMinkSumLSPlg((LineString)src, (Polygon)ref, doNeg, isRefConvex);
			}else if(src instanceof MultiPolygon) {
				return compMinkSumMultiPlgPlg((MultiPolygon)src, (Polygon)ref, doNeg, isRefConvex);
			}else if(src instanceof MultiLineString) {
				return compMinkSumMultiLSPlg((MultiLineString)src, (Polygon)ref, doNeg, isRefConvex);
			}else if(src instanceof GeometryCollection) {
				return compMinkSumGeometryCollection((GeometryCollection) src, ref, doNeg, isRefConvex);
			}else if(src instanceof Point) {
				return compMinkSumPoint((Point) src, ref, doNeg, isRefConvex);
			}
		}else if(ref instanceof LineString || ref instanceof LinearRing) {
			if(src instanceof Polygon) {
				return compMinkSumPlgLS((Polygon)src, (LineString)ref, doNeg);
			}else if(src instanceof LineString) {
				return compMinkSumLSLS((LineString)src, (LineString)ref, doNeg);
			}else if(src instanceof MultiPolygon) {
				return compMinkSumMultiPlgLS((MultiPolygon)src, (LineString)ref, doNeg);
			}else if(src instanceof MultiLineString) {
				return compMinkSumMultiLSLS((MultiLineString)src, (LineString)ref, doNeg);
			}else if(src instanceof GeometryCollection) {
				return compMinkSumGeometryCollection((GeometryCollection) src, ref, doNeg, isRefConvex);
			}else if(src instanceof Point) {
				return compMinkSumPoint((Point) src, ref, doNeg, isRefConvex);
			}
		}
		return src.getFactory().createPolygon();
	}
    //
	//trying out simulated multi-dispatch. Pity Java does not support it internally:( 
	//
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @param isRefConvex
	 * @return minkowski sum of src and ref (or an empty geometry)
	 */
	public static Geometry minkSum(Geometry src, Geometry ref, Boolean doNeg, Boolean isRefConvex) {
		if(src==null || ref == null) {
			return gf.createPolygon(); 
		}
		try {
			Class[] clzs = new Class[] {
					src.getClass(), ref.getClass(), Boolean.class, Boolean.class
			};
			Object[] objs = new Object[] {
					src, ref, doNeg, isRefConvex
			};
			Class clz = Class.forName("uk.osgb.algorithm.minkowski_sum.Minkowski_Sum"); 
			Method m = clz.getDeclaredMethod("minkSum", clzs);

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
		return src.getFactory().createPolygon(); // return an empty geometry
	}
	// wrappers for testing
	private static Geometry minkSum(Point src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumPoint(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(Point src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumPoint(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(Point src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumPoint(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(Polygon src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumPlgPlg(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(LineString src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumLSPlg(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(MultiPolygon src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumMultiPlgPlg(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(MultiLineString src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumMultiLSPlg(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(Polygon src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumPlgLS(src, ref, doNeg);
	}
	private static Geometry minkSum(Polygon src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumPlgLS(src, ref, doNeg);
	}
	private static Geometry minkSum(LineString src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumLSLS(src, ref, doNeg);
	}
	private static Geometry minkSum(LineString src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumLSLS(src, ref, doNeg);
	}
	private static Geometry minkSum(MultiPolygon src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumMultiPlgLS(src, ref, doNeg);
	}
	private static Geometry minkSum(MultiPolygon src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumMultiPlgLS(src, ref, doNeg);
	}
	private static Geometry minkSum(MultiLineString src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumMultiLSLS(src, ref, doNeg);
	}
	private static Geometry minkSum(MultiLineString src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumMultiLSLS(src, ref, doNeg);
	}
	private static Geometry minkSum(GeometryCollection src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumGeometryCollection(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(GeometryCollection src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumGeometryCollection(src, ref, doNeg, isRefConvex);
	}
	private static Geometry minkSum(GeometryCollection src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkSumGeometryCollection(src, ref, doNeg, isRefConvex);
	}
	//
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumPoint(Point src, Geometry ref, boolean doNeg, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			if(ref instanceof Polygon) {
				ref = negationPlg((Polygon) ref);
			}else if(ref instanceof LineString || ref instanceof LinearRing) {
				ref = negationLS((LineString) ref);
			}
		}
		AffineTransformation af = new AffineTransformation();
		af.translate(src.getX(), src.getY());
		return af.transform(ref);
	}
	
	public static Geometry compMinkSumMultiPoint(MultiPoint src, Geometry ref, boolean doNeg, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			if(ref instanceof Polygon) {
				ref = negationPlg((Polygon) ref);
			}else if(ref instanceof LineString || ref instanceof LinearRing) {
				ref = negationLS((LineString) ref);
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
	 * @param doNeg If negation is performed first on reference polygon (e.g. for collision detection purpose).
	 * For repeated computation, negation may be performed first before calling this method.
	 * @param isRefConvex whether the reference polygon is convex (if unknown, please use false)
	 * @return
	 */
	public static Geometry compMinkSumPlgPlg(Polygon src, Polygon ref, boolean doNeg, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkSumPlg(src, refCoords, isRefConvex);
	}
	//
	/** Minkowski sum of a multipolygon and a polygon
	 * @param src source multipolygon, may contain holes
	 * @param ref reference polygon, holes are ignored
	 * @param doNeg if negation is performed first on reference polygon (for repeated computation, negation may be performed first before calling this method)
	 * @param isRefConvex whether the reference polygon is convex (if unknown, please use false)
	 * @return
	 */
	public static Geometry compMinkSumMultiPlgPlg(MultiPolygon src, Polygon ref, boolean doNeg, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationPlg(ref);
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
	 * @param doNeg if negation is performed first on reference linestring (e.g. for collision detection purpose).
	 * For repeated computation, negation may be performed first before calling this method.
	 * @return
	 */
	public static Geometry compMinkSumPlgLS(Polygon src, LineString ref, boolean doNeg) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationLS(ref);
		}
		//
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkSumPlg(src, refCoords, false); // always treat LS ref as non-convex
	}
	//
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @return
	 */
	public static Geometry compMinkSumMultiPlgLS(MultiPolygon src, LineString ref, boolean doNeg) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationLS(ref);
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
	 * @param doNeg
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumLSPlg(LineString src, Polygon ref, boolean doNeg, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkSumLS(src, refCoords, isRefConvex);
	}
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @return
	 */
	public static Geometry compMinkSumLSLS(LineString src, LineString ref, boolean doNeg) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkSumLS(src, refCoords, false);

	}//
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @return
	 */
	public static Geometry compMinkSumMultiLSLS(MultiLineString src, LineString ref, boolean doNeg) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationLS(ref);
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
	 * @param doNeg
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumMultiLSPlg(MultiLineString src, Polygon ref, boolean doNeg, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationPlg(ref);
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
	 * @param doNeg
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkSumGeometryCollection(GeometryCollection src, Geometry ref, boolean doNeg, boolean isRefConvex) {
		if(src.isEmpty() || ref.isEmpty()) {
			return src.getFactory().createPolygon();
		}
		Geometry rlt = src.getFactory().createPolygon();
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Geometry part =  src.getGeometryN(i);
			Geometry sum1 = minkSum(part, ref, doNeg, isRefConvex);
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
	//
	/** Minkowski sum between a polygon (may contain holes) and a reference array of coordinates
	 * @param src Source polygon
	 * @param refCoords Coordinate sequence, which forms either an open LineString or a closed LinearRing 
	 * @param isRefConvex If refCoords is a closed LinearRing, whether it is convex or not
	 * @return
	 */
	private static Geometry compMinkSumPlg_backup(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		boolean isRefRing = false;
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring
			isRefRing = true;
		}
		Geometry sum = null;
		LineString extRing = src.getExteriorRing();
		Coordinate[] extCoords = extRing.getCoordinates();
		// expanding the shell of original geometry (on two directions)
		if(isRefRing) {
			sum = coordArrayVectorAddition(extCoords, refCoords, isRefConvex, src.getFactory());
		}else {
			sum = coordArrayVectorAddition(extCoords, refCoords, src.getFactory());
		}
		// if holes are outside original geometry, they should be kept 
		LinearRing[] shellHoles = getHoles((Polygon)sum);
		Vector<LinearRing> ringVec = null;
		
		if(shellHoles!=null) {
			ringVec = new Vector<LinearRing>(shellHoles.length);
			Polygon geomNoHole = gf.createPolygon(extCoords);// exterior boundary of result
			for(int i = 0; i < shellHoles.length; ++i) {
				Point ep = shellHoles[i].getEndPoint(); // a point on the hole boundary
				if(!geomNoHole.contains(ep)){// holes outside original geometry, generated by expansion, should be part of the result as a hole
					ringVec.add(shellHoles[i]);
				}
			}
			if(ringVec.size() == 0) {
				shellHoles = null;
			}else {
				shellHoles = new LinearRing[ringVec.size()];
				ringVec.toArray(shellHoles);
			}
		}else {
			ringVec = new Vector<LinearRing>();
		}
		// buffering hole bounaries in original geometry, if any and take interior boundaries
		Geometry sumHoles = null; // holes in the result, if any
		int numHoles = src.getNumInteriorRing();		
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				LineString holeRing = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeRing.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords);
				//Geometry sumHole = null;
				Polygon sumHole = null;
				
				if(isRefRing) {
					sumHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, isRefConvex, src.getFactory());
				}else {
					sumHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, src.getFactory());
				}
				if(sumHole!=null) {
					if(sumHole.getNumInteriorRing() > 0) {
						for(int i = 0; i < sumHole.getNumInteriorRing(); ++i){
							LineString ring = sumHole.getInteriorRingN(i);
							Point ep = ring.getEndPoint();
							if(holePlg.contains(ep)) {// not newly generated hole
								Polygon plgFromRing = gf.createPolygon(ring.getCoordinates());
								if(sumHoles==null) {
									sumHoles = plgFromRing;
								}else {
									try {
										sumHoles = sumHoles.union(plgFromRing);
									}catch(Exception e) {
										e.printStackTrace();
									}
								}
							}
						}
					}
				}
			}
		}
 
		LinearRing shell = getShellRing((Polygon)sum);
		if(sumHoles!=null) {
			// holes of sumHoles
			LinearRing[] holes = null;//
			//if(sumHoles.getGeometryType().equalsIgnoreCase("Polygon")) {
			if(sumHoles instanceof Polygon) {
				holes = getHoles((Polygon)sumHoles);
			}else if(sumHoles instanceof MultiPolygon) {
				holes = getHoles((MultiPolygon)sumHoles);
			}
			if(shellHoles!=null) {// merge holes generated by expansion to internal holes
				if(holes!=null) {
					for(int i = 0; i < holes.length; ++i) {
						ringVec.add(holes[i]);
					}
				}
				holes = new LinearRing[ringVec.size()];
				ringVec.toArray(holes);
			}
			return gf.createPolygon(shell, holes);
		}else {// no holes as result of expansion from original holes
			if(shellHoles==null) {// no holes
				return gf.createPolygon(shell.getCoordinates());
			}else {
				return gf.createPolygon(shell, shellHoles);
			}
		}		
	}
	//
	private static Geometry compMinkSumPlg_backup2(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
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
			Polygon geomNoHole = gf.createPolygon(extCoords);// exterior boundary of result
			for(int i = 0; i < shellHoles.length; ++i) {
				Point ep = shellHoles[i].getEndPoint(); // a point on the hole boundary
				if(!geomNoHole.contains(ep)){// holes outside original geometry, generated by expansion, should be part of the result as a hole
					if(ringVec == null) {
						ringVec = new Vector<LinearRing>();
					}
					ringVec.add(shellHoles[i]);
				}
			}
		}
		// buffering hole bounaries in original geometry, if any and take interior boundaries
		Geometry sumHoles = null; // holes in the result, if any
		int numHoles = src.getNumInteriorRing();		
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				LineString holeRing = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeRing.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords);
				Polygon sumHole = null;
				if(isRefRing) {
					sumHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, isRefConvex, src.getFactory());
				}else {
					sumHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, src.getFactory());
				}
				if(sumHole!=null) {
					if(sumHole.getNumInteriorRing() > 0) {
						for(int i = 0; i < sumHole.getNumInteriorRing(); ++i){
							LineString ring = sumHole.getInteriorRingN(i);
							Point ep = ring.getEndPoint();
							if(holePlg.contains(ep)) {// inside original hole, not newly generated hole
								if(ringVec==null) {
									ringVec = new Vector<LinearRing>();
								}
								ringVec.add(gf.createLinearRing(ring.getCoordinates()));
							}
						}
					}
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
		Geometry extSum = compMinkSumPlgNoHole(gf.createPolygon(extCoords), refCoords, isRefConvex);
		int numHoles = src.getNumInteriorRing();		
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				LineString holeRing = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeRing.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords);
				Geometry intDiff = compMinkDiffPlgNoHole(holePlg, refCoords, isRefConvex);
				extSum = extSum.difference(intDiff);
			}
		}
		return extSum;
	}

	private static Geometry compMinkSumPlgNoHole(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
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
				Point ep = shellHoles[i].getEndPoint(); // a point on the hole boundary
				if(!geomNoHole.contains(ep)){// holes outside original geometry, generated by expansion, should be part of the result as a hole
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
	 * 
	 *****************************************************************/
    //
	//trying out simulated multi-dispatch. Pity Java does not support it internally:( still has to support all subclasses manually
	// if src is not empty and ref is emtpy, difference is src; if src is empty and ref isnot empty, difference is empty set
	// 
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry minkDiff(Geometry src, Geometry ref, Boolean doNeg, Boolean isRefConvex) {
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
					src, ref, doNeg, isRefConvex
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
	private static Geometry minkDiff(Point src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return src.getFactory().createPolygon();	
	}
	private static Geometry minkDiff(MultiPoint src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return src.getFactory().createPolygon();	
	}
	private static Geometry minkDiff(Polygon src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffPlgPlg(src, ref, doNeg, isRefConvex);	
	}
	private static Geometry minkDiff(MultiPolygon src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffMultiPlgPlg(src, ref, doNeg, isRefConvex);	
	}
	private static Geometry minkDiff(Polygon src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffPlgLS(src, ref, doNeg);
	}
	private static Geometry minkDiff(Polygon src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffPlgLS(src, ref, doNeg);
	}
	private static Geometry minkDiff(MultiPolygon src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffMultiPlgLS(src, ref, doNeg);
	}
	private static Geometry minkDiff(MultiPolygon src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffMultiPlgLS(src, ref, doNeg);
	}
	private static Geometry minkDiff(GeometryCollection src, Polygon ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffGeometryCollection(src, ref, doNeg);
	}
	private static Geometry minkDiff(GeometryCollection src, LineString ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffGeometryCollection(src, ref, doNeg);
	}
	private static Geometry minkDiff(GeometryCollection src, LinearRing ref, Boolean doNeg, Boolean isRefConvex) {
		return compMinkDiffGeometryCollection(src, ref, doNeg);
	}
	/** Minkowski difference of two geometries. 
	 * @param src Source geometry, may be polygon or multipolygon
	 * @param ref Reference geometry, may be polygon or linestring/linearring
	 * @param doNeg
	 * @param isRefConvex
	 * @return null (not implemented yet)
	 */
	public static Geometry compMinkDiff(Geometry src, Geometry ref, boolean doNeg, boolean isRefConvex) {
		if(src==null) {
			return gf.createPolygon();
		}else if(ref==null) {
			return src;
		}
		if(ref instanceof Polygon) {
			if(src instanceof Polygon) {
				return compMinkDiffPlgPlg((Polygon)src, (Polygon)ref, doNeg, isRefConvex);
			}else if(src instanceof MultiPolygon) {
				return compMinkDiffMultiPlgPlg((MultiPolygon)src, (Polygon)ref, doNeg, isRefConvex);
			}else if(src instanceof GeometryCollection){
				return compMinkDiffGeometryCollection((GeometryCollection)src, ref, doNeg);
			}else {
				return src.getFactory().createPolygon();
			}
		}else if(ref instanceof LineString || ref instanceof LinearRing) {
			if(src instanceof Polygon) {
				return compMinkDiffPlgLS((Polygon)src, (LineString)ref, doNeg);
			}else if(src instanceof MultiPolygon) {
				return compMinkDiffMultiPlgLS((MultiPolygon)src, (LineString)ref, doNeg);
			}else if(src instanceof GeometryCollection){
				return compMinkDiffGeometryCollection((GeometryCollection)src, ref, doNeg);
			}else {
				return src.getFactory().createPolygon(); // unsupported source type
			}
		}else {
			return src.getFactory().createPolygon();
		}
	}
	//
	/**
	 * @param src Source polygon
	 * @param refCoords Coordinate sequence from reference geometry, negation if required is performed outside this method
	 * @param isRefConvex If the reference is a linear ring, whether it is convex
	 * @return
	 */
	private static Geometry compMinkDiffPlg_backup(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		boolean isRefRing = false;
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring so isRefConvex is considered
			isRefRing = true;
		}
		// bufferring exterior ring of source geometry
		Geometry sum = null;
		LineString extRing = src.getExteriorRing();
		Coordinate[] extCoords = extRing.getCoordinates();
		if(isRefRing) {
			sum = coordArrayVectorAddition(extCoords, refCoords, isRefConvex, src.getFactory());
		}else {
			sum = coordArrayVectorAddition(extCoords, refCoords, src.getFactory());
		}
		// 
		LinearRing[] shellHoles = getHoles((Polygon)sum);
		if(shellHoles==null) {
			return null;
		}
		// need to check if the holes are inside the outer boundary of the original geometry
		//Vector<LinearRing> ringVec = new Vector<LinearRing>(shellHoles.length);
		Vector<LinearRing> ringVec = new Vector<LinearRing>(shellHoles.length); // exterior boundary rings of mink diff
		
		Polygon geomNoHole = gf.createPolygon(src.getExteriorRing().getCoordinates());
		
		for(int i = 0; i < shellHoles.length; ++i) {
			Point ep = shellHoles[i].getEndPoint();
			if(geomNoHole.contains(ep)){ // the hole is inside the exterior ring of source geometry, not a hole generated by outwards expansion
				ringVec.add(shellHoles[i]);
			}
		}
		// buffering holes in source geometry
		int numHoles = src.getNumInteriorRing();
		// bueffering hole boundaries and take the exterior boundaries
		Geometry sumHoles = null; 
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				LineString holeRing = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeRing.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords); // hole polygon
				//Geometry sumHole = coordArrayVectorAddition(holeCoords, refCoords);
				//Geometry sumHole = null;
				Polygon sumHole = null;
				if(isRefRing) {
					sumHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, true, src.getFactory());
				}else {
					sumHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, src.getFactory());
				}
				if(sumHole!=null) {
					// exterior of sumHole
					LineString holeExtRing = sumHole.getExteriorRing();
					Polygon holeGeom = gf.createPolygon(holeExtRing.getCoordinates());
					if(sumHoles == null) {
						sumHoles = holeGeom;
					}else {
						try {
							sumHoles = sumHoles.union(holeGeom);
						}catch(Exception e) {
							e.printStackTrace();
						}
					}
					// check if there are newly generated hole which is part of the Mink difference
					if(sumHole.getNumInteriorRing() > 0) {
						for(int i = 0; i < sumHole.getNumInteriorRing(); ++i){
							LineString ring = sumHole.getInteriorRingN(i);
							Point ep = ring.getEndPoint();
							if(!holePlg.contains(ep)) { // a newly generated hole by boundary expansion, outside the source hole
								Coordinate[] ringCoords = ring.getCoordinates();
								if(ringCoords.length > 3 && ringCoords[0].equals2D(ringCoords[ringCoords.length-1])) {
									ringVec.add((LinearRing) ring);
								}// omit point
							}
						}
					}
				}
			}
		}
		//
		if(ringVec.size() == 0) {
			return null;
		}else {
			shellHoles = new LinearRing[ringVec.size()];
			ringVec.toArray(shellHoles);
		}
		//
		// there might be an issue here? one polygon (holes not added yet) contains another?
		//
		Geometry rltShell = null;
		if(shellHoles!=null && shellHoles.length > 1) {// multipolygon
			Polygon[] plgs = new Polygon[shellHoles.length];
			for(int i = 0; i < plgs.length; ++i) {
				plgs[i] = gf.createPolygon(shellHoles[i]);
			}
			rltShell = gf.createMultiPolygon(plgs);
		}else {// polygon
			rltShell = gf.createPolygon(shellHoles[0]);
		}
		// holes
		if(sumHoles!=null) {
			Geometry rltHoles = null;
			// shells of sumHoles
			LinearRing[] holeShells = null;//
			//if(sumHoles.getGeometryType().equalsIgnoreCase("Polygon")) {
			if(sumHoles instanceof Polygon) {
				holeShells = new LinearRing[1];
				holeShells[0] = getShellRing((Polygon)sumHoles);
				rltHoles = gf.createPolygon(holeShells[0]);
			//}else if(sumHoles.getGeometryType().equalsIgnoreCase("MultiPolygon")) {
			}else if(sumHoles instanceof MultiPolygon) {
				holeShells = getShellRings((MultiPolygon)sumHoles);
				if(holeShells!=null) {
					Polygon[] plgs = new Polygon[holeShells.length];
					for(int i = 0; i < plgs.length; ++i) {
						plgs[i] = gf.createPolygon(holeShells[i]);
					}
					rltHoles = gf.createMultiPolygon(plgs);
				}
			}
			if(rltHoles!=null) {
				return rltShell.difference(rltHoles);
			}else {
				return rltShell;
			}
		}else {// no holes in Mink diff
			return rltShell;
		}		
	}
	//
	/**
	 * @param src
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkDiffPlg_backup2(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
		boolean isRefRing = false;
		if(refCoords[0].equals2D(refCoords[refCoords.length-1])) { // reference is a linear ring so isRefConvex is considered
			isRefRing = true;
		}
	// bufferring exterior ring of source geometry
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
		if(shellHoles==null) {
			return gf.createPolygon(null, null);
		}
		//
	// need to check if the holes are inside the shell of the original geometry
	// holes could be generated from outwards expansion as well
		// external boundary of Mink Diff
		Vector<LinearRing> ringVec = new Vector<LinearRing>(shellHoles.length); 
		Polygon geomNoHole = gf.createPolygon(extCoords);
		for(int i = 0; i < shellHoles.length; ++i) {
			Point ep = shellHoles[i].getEndPoint();
			if(geomNoHole.contains(ep)){ // the hole is inside the exterior ring of source geometry, not a hole generated by outwards expansion
				ringVec.add(shellHoles[i]);
			}
		}
	// buffering holes in source geometry
		//  for each hole in source, buffer it and get exteriror as shell, and holes outside original hole as hole to assemble a geometry
		// union all hole geometry, then shell.difference(hole)
		int numHoles = src.getNumInteriorRing();
		// bueffering hole boundaries and take the exterior boundaries
		Geometry diffHoles = null; // union of all holes of the mink difference
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				
				LineString holeSrc = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeSrc.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords); // hole polygon
				Polygon diffHole = null;
				if(isRefRing) {
					diffHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, true, src.getFactory());
				}else {
					diffHole = (Polygon) coordArrayVectorAddition(holeCoords, refCoords, src.getFactory());
				}
				
				if(diffHole!=null) {
					// exterior of diffHole
					LineString holeExtRing = diffHole.getExteriorRing();
					Polygon holeGeom = gf.createPolygon(holeExtRing.getCoordinates());
					// check if there are newly generated hole which is part of the Mink difference
					if(diffHole.getNumInteriorRing() > 0) {
						Vector<LinearRing> holeVec = new Vector<LinearRing>();
						for(int i = 0; i < diffHole.getNumInteriorRing(); ++i){
							LineString ring = diffHole.getInteriorRingN(i);
							Point ep = ring.getEndPoint();
							if(!holePlg.contains(ep)) { // a newly generated hole by boundary expansion, outside the source hole
								holeVec.add(gf.createLinearRing(ring.getCoordinates()));
							}
						}
						if(holeVec.size() > 0) {
							LinearRing[] holes = new LinearRing[holeVec.size()];
							holeVec.toArray(holes);
							holeGeom = gf.createPolygon(gf.createLinearRing(holeExtRing.getCoordinates()), holes);
						}else {
							holeGeom = gf.createPolygon(holeExtRing.getCoordinates());
						}
					}else {
						holeGeom = gf.createPolygon(holeExtRing.getCoordinates());
					}
					if(diffHoles == null) {
						diffHoles = holeGeom;
					}else {
						try {
							diffHoles = diffHoles.union(holeGeom);
						}catch(Exception e) {
							e.printStackTrace();
						}
					}
				}
			}
		}
		//
		if(ringVec.size() == 0) {
			return null;
		}else {
			shellHoles = new LinearRing[ringVec.size()];
			ringVec.toArray(shellHoles);
		}
		// shell of result 
		Geometry rltShell = null;
		if(shellHoles!=null && shellHoles.length > 1) {// multipolygon
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
			if(diffHoles!=null) {
				return rltShell.difference(diffHoles);
			}else {// no holes in Mink diff
				return rltShell;
			}		
		}else {
			return gf.createPolygon(null, null);
		}
	}
	
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
		Geometry extDiff = compMinkDiffPlgNoHole(gf.createPolygon(extRing.getCoordinates()), refCoords, isRefConvex);
		int numHoles = src.getNumInteriorRing();		
		if(numHoles>0) {
			for(int k = 0; k < numHoles; ++k) {
				LineString holeSrc = src.getInteriorRingN(k);
				Coordinate[] holeCoords = holeSrc.getCoordinates();
				Polygon holePlg = gf.createPolygon(holeCoords); // hole polygon
				Geometry intSum = compMinkSumPlgNoHole(holePlg, refCoords, isRefConvex);
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
	 * @param src polygon without hole
	 * @param refCoords
	 * @param isRefConvex
	 * @return
	 */
	private static Geometry compMinkDiffPlgNoHole(Polygon src, Coordinate[] refCoords, boolean isRefConvex) {
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
			Point ep = shellHoles[i].getEndPoint();
			if(geomNoHole.contains(ep)){ // the hole is inside the original source geometry, not a hole generated by outwards expansion
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
	 * @param doNeg whether ref should be negated (for repeated calls, ref may be negated before calling this method). If ref is symmetric to the coordinate origin, there is no need to perform negation either.
	 * @PARAM isRefConvex If you know the ref polygon is convex, set this to true will improve performace and robustness.
	 * @return
	 */
	public static Geometry compMinkDiffPlgPlg(Polygon src, Polygon ref, boolean doNeg, boolean isRefConvex) {
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationPlg(ref);
		}
		Coordinate[] refCoords = ref.getExteriorRing().getCoordinates();
		return compMinkDiffPlg(src, refCoords, isRefConvex);
	}
	
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @return
	 */
	public static Geometry compMinkDiffPlgLS(Polygon src, LineString ref, boolean doNeg) {
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationLS(ref);
		}
		Coordinate[] refCoords = ref.getCoordinates();
		return compMinkDiffPlg(src, refCoords, false);
	}
	/**
	 * @param src
	 * @param ref
	 * @param doNeg
	 * @param isRefConvex
	 * @return
	 */
	public static Geometry compMinkDiffMultiPlgPlg(MultiPolygon src, Polygon ref, boolean doNeg, boolean isRefConvex) {
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationPlg(ref);
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
	 * @param doNeg
	 * @return
	 */
	public static Geometry compMinkDiffMultiPlgLS(MultiPolygon src, LineString ref, boolean doNeg) {
		if(doNeg) {// for geometry symmetric on coordinate origin, the negatied geometry is the same as original
			ref = negationLS(ref);
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
	public static Geometry compMinkDiffGeometryCollection(GeometryCollection src, Geometry ref, boolean doNeg) {
		Geometry rlt = null;
		int numParts = src.getNumGeometries();
		for(int i = 0; i < numParts; ++i) {
			Geometry geom = src.getGeometryN(i);
			Geometry diff1 = minkDiff(geom, ref, doNeg, false);
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
	 * @return negated geometry (null if type is not supported yet)
	 */
	public static Geometry negation(Geometry geom) {
		/*
		String gType = geom.getGeometryType();
		if(gType.equalsIgnoreCase("POLYGON")) {
			return negationPlg((Polygon)geom);
		}else if(gType.equalsIgnoreCase("LINESTRING") ||gType.equalsIgnoreCase("LINEARRING")) {
			return negationLS((LineString)geom);
		}
		*/
		if(geom instanceof Polygon) {
			return negationPlg((Polygon)geom);
		}else if(geom instanceof LineString || geom instanceof LinearRing) {
			return negationLS((LineString)geom);
		}
		return null;
	}
	//
	/** negation of JTS polygon
	 * @param plg
	 * @return
	 */
	public static Polygon negationPlg(Polygon plg){
		Coordinate[] extNeg = negation(plg.getExteriorRing().getCoordinates());
		GeometryFactory gf = plg.getFactory();
		int numHoles = plg.getNumInteriorRing();
		if(numHoles > 0) {
			LinearRing[] holes = new LinearRing[numHoles];
			for(int i = 0; i < numHoles; ++i) {
				Coordinate[] holeNeg = negation(plg.getInteriorRingN(i).getCoordinates());
				holes[i] = gf.createLinearRing(holeNeg);
			}
			return gf.createPolygon(gf.createLinearRing(extNeg), holes);
		}else {
			return gf.createPolygon(extNeg);
		}
	}
	public static LineString negationLS(LineString ls){
		Coordinate[] coordsNeg = negation(ls.getCoordinates());
		return ls.getFactory().createLineString(coordsNeg);
	}
	/** negation of a Coordinate array
	 * @param coords
	 * @return
	 */
	public static Coordinate[] negation(Coordinate[] coords) {
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
	 * @param refCoords Coordinate array of the reference polygon without holes
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
	public static Geometry segLSAddition(Coordinate[] refCoords, Coordinate segSp, Coordinate segEp, GeometryFactory gf) {
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
	/** vector addition of a segment on reference linestring and a segment on the source geometry
	 *  Convexhull is used to construct result to improve robustness (avoiding self-intersection, hopefully).
	 * @param segSp
	 * @param segEp
	 * @param refSp
	 * @param refEp
	 * @return a polygon (empty polygon is the sum is 0D or 1D)
	 */
	private static Geometry segVectorAddition(Coordinate segSp, Coordinate segEp, Coordinate refSp, Coordinate refEp, GeometryFactory gf){
		Coordinate[] coords = new Coordinate[4];
		coords[0] = new Coordinate(refSp.x + segSp.x, refSp.y + segSp.y);
		coords[1] = new Coordinate(refEp.x + segSp.x, refEp.y + segSp.y);
		coords[2] = new Coordinate(refEp.x + segEp.x, refEp.y + segEp.y);
		coords[3] = new Coordinate(refSp.x + segEp.x, refSp.y + segEp.y);
		ConvexHull ch = new ConvexHull(coords, gf);
		Geometry hull = ch.getConvexHull();
		//if(hull.getNumPoints() > 2) {
		if(hull.getDimension() > 1) {
			return hull; // a polygon
		}else {
			//return null; // collinear, ignored
			return gf.createPolygon(); // return an empty polygon?
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
				Geometry sum1 = segLSAddition(refCoords, segSp, segEp, gf);
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
	
	public static void setGeometryFactory(GeometryFactory geomFactory) {
		gf = geomFactory;
	}
	
	public static void main(String[] args) {

/*		Coordinate[] refCoords01 = new Coordinate[9];
		refCoords01[0] = new Coordinate(0.5, 0);
		refCoords01[1] = new Coordinate(0.125, 0.125);
		refCoords01[2] = new Coordinate(0, 0.5);
		refCoords01[3] = new Coordinate(-0.125, 0.125);
		refCoords01[4] = new Coordinate(-0.5, 0);
		refCoords01[5] = new Coordinate(-0.125, -0.125);
		refCoords01[6] = new Coordinate(0, -0.5);
		refCoords01[7] = new Coordinate(0.125, -0.125);
		refCoords01[8] = new Coordinate(0.5, 0);
*/
		Coordinate[] refCoords00 = new Coordinate[5];
		refCoords00[0] = new Coordinate(0.5, -0.5);
		refCoords00[1] = new Coordinate(0.5, 0.5);
		refCoords00[2] = new Coordinate(-0.5, 0.5);
		refCoords00[3] = new Coordinate(-0.5, -0.5);
		refCoords00[4] = new Coordinate(0.5, -0.5);
		
		//Polygon ref01 = gf.createPolygon(refCoords01);
		Geometry ref00 = gf.createPolygon(refCoords00);
		AffineTransformation at = new AffineTransformation();
		at.scale(10, 10);
		//ref01 = (Polygon) at.transform(ref01);
		ref00 = at.transform(ref00);
/*		
		Coordinate[] refCoords02 = new Coordinate[4];
		refCoords02[0] = new Coordinate(-0.25, -0.25);
		refCoords02[1] = new Coordinate(0, -0.25);
		refCoords02[2] = new Coordinate(0, 0.25);
		refCoords02[3] = new Coordinate(0.25, 0.25);
		LineString ref02 = gf.createLineString(refCoords02);
		
		Coordinate[] refCoords03 = new Coordinate[5];
		refCoords03[0] = new Coordinate(-0.25, -0.25);
		refCoords03[1] = new Coordinate(0.25, -0.25);
		refCoords03[2] = new Coordinate(0.25, 0.25);
		refCoords03[3] = new Coordinate(-0.25, 0.25);
		refCoords03[4] = new Coordinate(-0.25, -0.25);
		//LineString ref03 = gf.createLineString(refCoords02);
		LinearRing ref03 = gf.createLinearRing(refCoords03);

		Coordinate[] plgCoords = new Coordinate[6];
		plgCoords[0] = new Coordinate(4, 2);
		plgCoords[1] = new Coordinate(14, 2);
		plgCoords[2] = new Coordinate(12, 6);
		plgCoords[3] = new Coordinate(9, 5);
		plgCoords[4] = new Coordinate(6, 6);
		plgCoords[5] = new Coordinate(4, 2);
		LinearRing extRing = gf.createLinearRing(plgCoords);
		
		LinearRing[] holes = new LinearRing[1];
		Coordinate[] holeCoords = new Coordinate[5];
		holeCoords[0] = new Coordinate(8.25, 3.25);
		holeCoords[1] = new Coordinate(9.75, 3.25);
		holeCoords[2] = new Coordinate(9.75, 4.75);
		holeCoords[3] = new Coordinate(8.25, 4.75);
		holeCoords[4] = new Coordinate(8.25, 3.25);
		
		holes[0] = gf.createLinearRing(holeCoords);
		//Polygon plg = gf.createPolygon(extRing, holes);
		Geometry plg = gf.createPolygon(extRing, holes);
		Coordinate[] lsCoords = new Coordinate[3];
		lsCoords[0] = new Coordinate(5, 5);
		lsCoords[1] = new Coordinate(10, 5);
		lsCoords[2] = new Coordinate(10, 15);
		LineString ls = gf.createLineString(lsCoords);
		
		//Geometry sum = compMinkSumPlgPlg(plg, ref01, false, false);
		Geometry sum = minkSum(plg, ref01, false, false);
		System.out.println("Minkowski sum of two polygons: " + sum.toText());
		//sum = compMinkSum(plg, ref02, false, false);
		sum = minkSum(plg, ref02, false, false);
		System.out.println("Minkowski sum of a polygon and a linestring: " + sum.toText());
		sum = minkSum(plg, ref03, false, false);
		System.out.println("Minkowski sum of a polygon and a linearring: " + sum.toText());

		//sum = compMinkSum(ls, ref01, false, false);
		sum = minkSum(ls, ref01, false, false);
		System.out.println("Minkowski sum a linestring and a polygon " + sum.toText());
		//sum = compMinkSum(ls, ref02, false, false);
		sum = minkSum(ls, ref02, false, false);
		System.out.println("Minkowski sum of two linestring: " + sum.toText());
		
		//Geometry diff = compMinkDiffPlgPlg(plg, ref01, false, false);
		Geometry diff = minkDiff(plg, ref01, false, false);
		System.out.println("Minkowski difference of two polygons: " + diff.toText());
		//diff = compMinkDiff(plg, ref02, false, false);
		diff = minkDiff(plg, ref02, false, false);
		System.out.println("Minkowski difference of a polygon and a linestring: " + diff.toText());
*/
		String wktStr1 = "POLYGON ((350740.0283040225 276585.561347613, 351092.4537635583 276589.7702354695, 351088.0990496299 276526.9379345023, 351042.06350238656 276526.9379345023, 351033.6651255246 276558.9761869757, 350982.34171136824 276550.2667591189, 350983.2748643529 276504.54226287047, 351035.2203804991 276503.2980588909, 351041.4414003968 276522.272169579, 351087.4769476401 276521.96111858415, 351086.23274366057 276479.3471322846, 350750.028484311 276476.2736630311, 350740.0283040225 276585.561347613), (350797.8864899776 276556.9894039315, 350809.31526745023 276512.70289122505, 350886.538004943 276512.6295887375, 350886.2269539481 276537.2026173336, 350908.0005235902 276537.2026173336, 350911.11103353905 276515.1179966966, 350937.5503681045 276514.18484371196, 350941.28298004315 276556.7988300115, 350908.62262558 276558.04303399107, 350908.93367657484 276541.55733126204, 350885.60485195834 276542.1794332518, 350882.8053930043 276560.5314419501, 350797.8864899776 276556.9894039315))";
		String wktStr2 = "POLYGON ((350746.87142591 276458.652541699, 350991.8758429792 276462.9383332512, 350989.5621953245 276414.2096860992, 351040.4724679055 276414.8582246034, 351044.019640198 276459.36684029107, 351099.734930377 276450.08095859457, 351093.30624304863 276354.3649472614, 351041.16244582983 276356.5078430375, 351043.30534160597 276408.65164025634, 350990.4472457951 276406.5087444802, 350988.30435001897 276350.0791557092, 350756.8716061985 276349.3648571171, 350746.87142591 276458.652541699), (350804.72961186507 276430.0805980175, 350816.1583893377 276385.79408531106, 350936.16055280017 276384.365488127, 350937.67911498697 276405.07069767715, 350922.58887955145 276410.0802374404, 350906.8743105266 276397.2228627837, 350889.73114431766 276407.22304307227, 350893.36231719883 276418.0414677615, 350924.2759858998 276415.663493246, 350939.0177471683 276418.65182054485, 350941.87494153646 276435.7949867538, 350804.72961186507 276430.0805980175))";
		WKTReader rd = new WKTReader();
		try {
			Geometry geom1 = rd.read(wktStr1);
			Geometry sum = minkSum(geom1, ref00, false, false);
			System.out.println(sum.toText());
			Geometry sumDiff = minkDiff(sum, ref00, false, false);
			Geometry geom2 = rd.read(wktStr2);
			Geometry diff = minkDiff(geom2, ref00, false, false);
			System.out.println(diff.toText());
			Geometry diffSum = minkSum(diff, ref00, false, false);
			String sumStr = sumDiff.toText();
			System.out.println(sumStr);
			String diffStr = diffSum.toText();
			System.out.println(diffStr); 
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
