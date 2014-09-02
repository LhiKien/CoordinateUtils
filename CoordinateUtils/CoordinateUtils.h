//
//  CoordinateUtils.h
//  CoordinateUtils
//
//  Created by damingdan on 14/8/27.
//  Copyright (c) 2014年 Kingoit. All rights reserved.
//

#import <Foundation/Foundation.h>

typedef struct _MapParameter {
    double deltaX;
    double deltaY;
    double deltaZ;
    double rotateX;
    double rotateY;
    double rotateZ;
    double k;
    double offsetX;
    double offsetY;
}MapParameter;

/**
 参心大地坐标系中的点
 */
typedef struct _GeodeticPoint {
    double latitude; // 纬度
    double longitude; // 经度
    double height; // 海拔
}GeodeticPoint;

/**
 空间直角坐标系中的点
 */
typedef struct _RectangularPlanePoint {
    double x;
    double y;
    double z;
}RectangularPlanePoint;

typedef struct _GaussPoint {
    double x;
    double y;
    double h;
}GaussPoint;

/**
 西安80坐标系定义的椭球参数：
 a ＝ 6378140
 b = 6356755.2882
 f = 1/298.257
 e^2 = 0.00669438499959
 
 WGS-84 椭球参数:
 a = 6378137
 b = 6356752.314
 f = 1/298.257223563
 e^2 = 0.006694379989
 
 克氏椭球参数(即北京54坐标系):
 a = 6378245
 b = 6356863.019
 f = 1/298.3
 e^2 = 0.006693421623
 
 注：a, b ,f, e分别表示椭球体的长半轴，短半轴，扁率和第一偏心率, 单位米
 */

@interface Ellipsoid : NSObject
/**
 椭球的长半轴
 */
@property(nonatomic) double a;

/**
 椭球的短半轴
 */
@property(nonatomic) double b;

/**
 扁率
 */
@property(nonatomic) double f;

/**
 第一偏心率
 */
@property(nonatomic) double e2;

/**
 第二偏心率
 */
@property(nonatomic) double ee;

+ (Ellipsoid*) WGS84;

+ (Ellipsoid*) Xian80;

+ (Ellipsoid*) BeiJin54;
@end


@interface CoordinateUtils : NSObject

/**
 参心大地坐标转换为参心空间直角坐标
 坐标转换的公式为：
 X = (N+H)*cosB*cosL
 Y = (N+H)*cosB*sinL
 Z = [N*(1-e^2) + H]*sinB
 X, Y, Z分别代表空间直角坐标系中的点的三个参数
 B, L, H分别表示纬度， 经度， 和海拔
 N表示椭球面卯酉圈的曲率半径，e^2为椭球的第一偏心率，a、b 椭球的长短半径，f 椭球扁率，W为第一辅助系数
 e^2 = (a*a - b*b)/(a*a) = 1 - pow(b/a, 2)
 
 W = sqrt(1 - e^2*sin(B)^2)
 
 N = a/W
 
 @param point 大地坐标点
 @param ellipsoid 坐标转换所在的椭球
 
 @return 相对应的空间直角坐标系中的点
 */
- (RectangularPlanePoint) rectangularPlaneCoordinateFromGeodeticCoordinates:(GeodeticPoint) point inEllipsoid:(Ellipsoid*) ellipsoid;


/**
 空间直角坐标系转大地坐标系
 L = arctan(Y/X)
 B的计算通常采用牛顿迭代法来近似计算，理论上可以达到任意想要的精度
 B(0) = arctan(Z/sqrt(X^2 + Y^2))
 B(i+1) = arctan((Z + N(i)*e^2*sin(B(i)))/sqrt(X^2 + Y^2))
 H = Z/sin(B) - N*(1 - e^2)
 
 @param point 空间直角坐标系中的点
 @param ellipsoid 计算所基于的椭球定义
 @return 相对应的大地坐标系中的点
 */
- (GeodeticPoint) geodeticCoordinatesFromRectangularPlaneCoordinate:(RectangularPlanePoint) point inEllipsoid:(Ellipsoid*) ellipsoid;

/**
 高斯投影正算，计算大地坐标系（就是经纬度坐标系）中的点投影到2D平面坐标系中的点
 
 测试数据：
 正常的经纬度值：  纬度：32  经度：121
 北京65投影结果：  X：310994 Y：3543664
 WGS84投影结果：  X：310997 Y：3543601
 
 @param geoPoint 大地坐标点，即经纬度坐标
 @param ellipsoid 计算所基于的椭球定义
 @param zoneWidth 带宽，一般为3度分带活着6度分带
 @return 高斯投影所得到的平面坐标点
 */
- (GaussPoint) gaussProjection:(GeodeticPoint) geoPoint inEllipsoid:(Ellipsoid*) ellipsoid withZoneWidth:(int) zoneWidth;


/**
 计算高斯投影反算中非常重要的参数Bf（即底点纬度)
 @param x 高斯投影坐标系中的x坐标
 @param ellipsoid 计算所基于的椭球定义
 @return 地点纬度的值
 */
- (double) calculateBfWithX:(double)x inEllipsoid:(Ellipsoid*)ellipsoid;


/**
 高斯投影反算, 带带号， 直接给出所处的带号
 @param gaussPoint  所需要计算转化的高斯平面投影点
 @param ellipsiod 计算所基于的椭球定义
 @param zoneWidth 带宽，一般是3度分带或者6度分带
 @return 经过高斯反算所得到的经纬度坐标 @c GeodeticPoint
 */
- (GeodeticPoint) gaussProjectionReverse:(GaussPoint) gaussPoint inEllipsoid:(Ellipsoid*) ellipsiod withZoneWidth:(int) zoneWidth;

/**
 高斯投影反算, 不带带号， 直接给出所处的带号
 @param gaussPoint  所需要计算转化的高斯平面投影点
 @param ellipsiod 计算所基于的椭球定义
 @param zoneWidth 带宽，一般是3度分带或者6度分带
 @param projectionNumber 带号， 因为坐标中没有给出带号，所以这里需要自己指定带号
 @return 经过高斯反算所得到的经纬度坐标 @c GeodeticPoint
 */
- (GeodeticPoint) gaussProjectionReverse:(GaussPoint) gaussPoint inEllipsoid:(Ellipsoid*) ellipsiod withZoneWidth:(int) zoneWidth inProjectNumber:(int) projectionNumber;

/**
 84坐标点转80坐标点
 */
- (GaussPoint) Xian80PointFromWGS84Point:(GeodeticPoint) point;
@end
