//
//  CoordinateUtils.m
//  CoordinateUtils
//
//  Created by damingdan on 14/8/27.
//  Copyright (c) 2014年 Kingoit. All rights reserved.
//

#import "CoordinateUtils.h"
#import <math.h>

@implementation CoordinateUtils


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
- (RectangularPlanePoint) rectangularPlaneCoordinateFromGeodeticCoordinates:(GeodeticPoint) point inEllipsoid:(Ellipsoid*) ellipsoid {
    RectangularPlanePoint recttangularPlanePoint;
    
    double latitude = [self translateAngle:point.latitude];
    double longitude = [self translateAngle:point.longitude];
    double height = point.height;
    
    double W = sqrt(1 - ellipsoid.e2 * pow(sin(latitude), 2));
    double N = ellipsoid.a/W;
    
    recttangularPlanePoint.x = (N + height)*cos(latitude)*cos(longitude);
    
    recttangularPlanePoint.y = (N + height)*cos(latitude)*sin(longitude);
    
    recttangularPlanePoint.z = (N*(1 - ellipsoid.e2) + height)*sin(latitude);
    
    return recttangularPlanePoint;
}

/**
 用布尔莎七参数公式计算两个不同的空间直角坐标系之间的转化
 */
- (RectangularPlanePoint) rectangularPlaneTanslate:(RectangularPlanePoint) point withParameter:(MapParameter) parameter{
    RectangularPlanePoint toPoint;
    
    double deltaX = parameter.deltaX;
    double deltaY = parameter.deltaY;
    double deltaZ = parameter.deltaZ;
    //标准参数中给出的单位是秒，这里要转成弧度
    double rotateX = [self translateAngle:parameter.rotateX/3600.0];
    double rotateY = [self translateAngle:parameter.rotateY/3600.0];
    double rotateZ = [self translateAngle:parameter.rotateZ/3600.0];
    // 标准参数中给出的值的单位是ppm，百万分之一米
    double k = 1.0 / parameter.k / 10000000.0;
    
    //布尔莎七参数公式的另一种形式
//    toPoint.x = deltaX + (1 + k)*(point.x + rotateZ*point.y - rotateY*point.z);
//    toPoint.y = deltaY + (1 + k)*(-rotateZ*point.x+ point.y - rotateX*point.z);
//    toPoint.z = deltaZ + (1 + k)*(point.x*rotateY - rotateX*point.y - point.z);
    
    toPoint.x = deltaX + (1 + k)*point.x + rotateZ*point.y - rotateY*point.z;
    toPoint.y = deltaY + (1 + k)*point.y - rotateZ*point.x + rotateX*point.z;
    toPoint.z = deltaZ + (1 + k)*point.z + rotateY*point.x - rotateX*point.y;
    
    return toPoint;
}

/**
 空间直角坐标系转大地坐标系
 L = arctan(Y/X)
 B的计算通常采用牛顿迭代法来近似计算，理论上可以达到任意想要的精度
 B(0) = arctan(Z/sqrt(X^2 + Y^2))
 B(i+1) = arctan((Z + N(i)*e^2*sin(B(i)))/sqrt(X^2 + Y^2))
 H = Z/sin(B) - N*(1 - e^2)
 */
- (GeodeticPoint) geodeticCoordinatesFromRectangularPlaneCoordinate:(RectangularPlanePoint) point inEllipsoid:(Ellipsoid*) ellipsoid {
    GeodeticPoint geodeticPoint;
    
    geodeticPoint.longitude = atan(point.y/point.x);
    if(geodeticPoint.longitude < 0) {
        geodeticPoint.longitude += M_PI;
    }
    
    double b0 = atan(point.z/sqrt(point.x*point.x + point.y*point.y));// 纬度迭代的起始值
    double delta = M_PI / (180000 * 3600);// 迭代的精度控制，这个值越小，精度越高
    double N = ellipsoid.a/sqrt(1 - ellipsoid.e2 * pow(sin(b0), 2));
    double b1 = 0;
    
    while (b0 - b1 >= delta) {
        b1 = b0;
        b0 = atan((point.z + N*ellipsoid.e2*sin(b1))/sqrt(point.x*point.x + point.y*point.y));
        N = ellipsoid.a/sqrt(1 - ellipsoid.e2 * pow(sin(b0), 2));
    }
    geodeticPoint.latitude = b0;
    if (geodeticPoint.latitude < 0) {
        geodeticPoint.latitude += M_PI;
    }
    
    geodeticPoint.height = point.z/sin(geodeticPoint.latitude) - N*(1 - ellipsoid.e2);
    // 将他们的弧度转成角度
    geodeticPoint.longitude = [self tanslateDegree:geodeticPoint.longitude];
    geodeticPoint.latitude = [self tanslateDegree:geodeticPoint.latitude];
    return geodeticPoint;
}

/**
 弧度转角度
 */
- (double) tanslateDegree:(double) degree {
    return degree*180/M_PI;
}

/**
 角度转弧度
 */
- (double) translateAngle:(double) angle {
    return angle/180.0*M_PI;
}

/**
 将度分秒转化成度数
 */
- (double) translateDMSToAngle:(double) dms{
    int Deg,Min;
    double Sec;
    Deg=(int)dms;
    Min=(int)((dms-Deg)*100);
    Sec=((dms-Deg)*100-Min)*100;
    return (Deg+Min/60.0+Sec/3600.0);
}


/**
 高斯投影正算，计算大地坐标系（就是经纬度坐标系）中的点投影到2D平面坐标系中的点
 
 测试数据：
 正常的经纬度值：  纬度：32  经度：121
 北京65投影结果：  X：310994 Y：3543664
 WGS84投影结果：  X：310997 Y：3543601
 */
- (GaussPoint) gaussProjection:(GeodeticPoint) geoPoint inEllipsoid:(Ellipsoid*) ellipsoid withZoneWidth:(int) zoneWidth{
    GaussPoint point;
    
    // 计算带号
    int projectionNumber;
    // 计算代号中央经线
    double longtitude0;
    if (zoneWidth == 3) {
        projectionNumber = (int)(geoPoint.longitude/zoneWidth);
        longtitude0 = projectionNumber*zoneWidth;
    }else if(zoneWidth == 6) {
        projectionNumber = (int)(geoPoint.longitude/zoneWidth) + 1;
        longtitude0 = projectionNumber*zoneWidth - zoneWidth/2.0;
    }
    
    
    // 将角度转化为弧度
    longtitude0 = [self translateAngle:longtitude0];
    double longtitude = [self translateAngle:geoPoint.longitude];
    double lantitude = [self translateAngle:geoPoint.latitude];
    
    // 计算椭球面卯酉圈的曲率半径
    double N = ellipsoid.a/sqrt(1 - ellipsoid.e2 * pow(sin(lantitude), 2));
    
    // 计算公式中的参数T
    double T = pow(tan(lantitude), 2);
    
    // 计算公式中用到的参数C
    double C = ellipsoid.ee * pow(cos(lantitude), 2);
    
    // 计算公式中用到的参数A
    double A = (longtitude - longtitude0) * cos(lantitude);
    
    // 算克拉索夫斯基椭球中子午弧长M，数学中常用X表示
    double M = ellipsoid.a * ((1 - ellipsoid.e2/4 - 3*pow(ellipsoid.e2, 2)/64 - 5*pow(ellipsoid.e2, 3)/256)*lantitude - (3*ellipsoid.e2/8 + 3*pow(ellipsoid.e2, 2)/32 + 45*pow(ellipsoid.e2, 3)/1024)*sin(2*lantitude) + (15*pow(ellipsoid.e2, 2)/256 + 45*pow(ellipsoid.e2, 3)/1024)*sin(4*lantitude) - (35*pow(ellipsoid.e2, 3)/3072)*sin(6*lantitude));
    
    //因为是以赤道为Y轴的，与我们南北为Y轴是相反的，所以xy与高斯投影的标准xy正好相反;
    point.x = N*(A + (1 - T + C)*pow(A, 3)/6 + (5 - 18*T + pow(T, 2) + 14*C - 58*ellipsoid.ee)*pow(A, 5)/120);
    
    point.y = M + N*tan(lantitude)*(A*A/2 + (5 - T + 9*C + 4*C*C)*pow(A, 4)/24 + (61 - 58*T + T*T + 270*C - 330*ellipsoid.ee)*pow(A, 6)/720);
    
    // 在我国 坐标都是正的， 坐标的最大值（在赤道上）约为330km。为了避免出现负的横坐标，可在横坐标上加上500 000m。
    point.x += projectionNumber*1000000L + 500000L;
    point.h = geoPoint.height;
    
    return point;
}

/**
 计算高斯投影反算中非常重要的参数Bf（即底点纬度）
 */
- (double) calculateBfWithX:(double)x inEllipsoid:(Ellipsoid*)ellipsoid {
    double e2 = ellipsoid.e2;
    double e4 = pow(ellipsoid.e2, 2);
    double e6 = pow(ellipsoid.e2, 3);
    double e8 = pow(ellipsoid.e2, 4);
    double e10 = pow(ellipsoid.e2, 5);
    double e12 = pow(ellipsoid.e2, 6);
    double e14 = pow(ellipsoid.e2, 7);
    double e16 = pow(ellipsoid.e2, 8);
    
    double c0 = 1 + e2 / 4 + 7 * e4 / 64 + 15 * e6 / 256 + 579 * e8 / 16384 + 1515 * e10 / 65536 + 16837 * e12 / 1048576 + 48997 * e14 / 4194304 + 9467419 * e16 / 1073741824;
    
    c0 = ellipsoid.a/c0;
    
    double b0 = x/c0;
    
    double d1 = 3 * e2 / 8 + 45 * e4 / 128 + 175 * e6 / 512 + 11025 * e8 / 32768 + 43659 * e10 / 131072 + 693693 * e12 / 2097152 + 10863435 * e14 / 33554432;
    
    double d2 = -21 * e4 / 64 - 277 * e6 / 384 - 19413 * e8 / 16384 - 56331 * e10 / 32768 - 2436477 * e12 / 1048576 - 196473 * e14 / 65536;
    
    double d3 = 151 * e6 / 384 + 5707 * e8 / 4096 + 53189 * e10 / 163840 + 4599609 * e12 / 655360 + 15842375 * e14 / 1048576;
    
    double d4 = -1097 * e8 / 2048 - 1687 * e10 / 640 - 3650333 * e12 / 327680 - 114459079 * e14 / 27525120;
    
    double d5 = 8011 * e10 / 1024 + 874457 * e12 / 98304 + 216344925 * e14 / 3670016;
    
    double d6 = -682193 * e12 / 245760 - 46492223 * e14 / 1146880;
    
    double d7 = 36941521 * e14 / 3440640;
    
    double bf;
    bf=b0+sin(2*b0)*(d1+sin(b0)*sin(b0)*(d2+sin(b0)*sin(b0)*(d3+sin(b0)*sin(b0)*(d4+sin(b0)*sin(b0)*(d5+sin(b0)*sin(b0)*(d6+d7*sin(b0)*sin(b0)))))));
    return bf;
}

/**
 高斯投影反算, 计算带带号的坐标
 */
- (GeodeticPoint) gaussProjectionReverse:(GaussPoint) gaussPoint inEllipsoid:(Ellipsoid*) ellipsiod withZoneWidth:(int) zoneWidth {
    // 计算带号
    int projectionNumber = (int)(gaussPoint.x/1000000L);
    gaussPoint.x -= projectionNumber*1000000L;
    return [self gaussProjectionReverse:gaussPoint inEllipsoid:ellipsiod withZoneWidth:zoneWidth inProjectNumber:projectionNumber];
}

/**
 高斯投影反算, 不带带号， 直接给出所处的带号
 */
- (GeodeticPoint) gaussProjectionReverse:(GaussPoint) gaussPoint inEllipsoid:(Ellipsoid*) ellipsiod withZoneWidth:(int) zoneWidth inProjectNumber:(int) projectionNumber{
    GeodeticPoint point;
    // 计算代号中央经线
    double longtitude0;
    if(zoneWidth == 3) {
        longtitude0 = projectionNumber*zoneWidth;
    }else if(zoneWidth == 6) {
        longtitude0 = (projectionNumber - 1)*zoneWidth + 0.5*zoneWidth;
    }
    
    longtitude0 = [self translateAngle:longtitude0];
    
    // 去除代号和偏移值
    double x = gaussPoint.y;
    double y = gaussPoint.x - 500000L;
    
    // 计算反算公式中的各个符号的值
    double Bf = [self calculateBfWithX:x inEllipsoid:ellipsiod];
    double Tf = tan(Bf);
    double Nf2 = ellipsiod.ee*cos(Bf)*cos(Bf);
    double Wf = sqrt(1 - ellipsiod.e2*sin(Bf)*sin(Bf));
    double Mf = ellipsiod.a*(1 - ellipsiod.e2)/pow(Wf, 3);
    double Nf = ellipsiod.a/Wf;
    
    point.longitude = longtitude0 +  y/(Nf* cos(Bf)) - (1 + 2*Tf*Tf + Nf2)*pow(y, 3)/(6*pow(Nf, 3)*cos(Bf)) + (5 + 28*Tf*Tf + 24*pow(Tf, 4))*pow(5, 5)/(120*pow(Nf, 5)*cos(Bf));
    
    point.latitude = Bf - Tf*y*y/(2*Mf*Nf) + Tf*(5 + 3*Tf + Nf2 - 9*Nf2*Tf*Tf)*pow(y, 4)/(24*Mf*pow(Nf, 3)) - Tf*(61 + 90*Tf*Tf + 45*pow(Tf, 4))*pow(y, 6)/(720*Mf*pow(Nf, 5));
    
    point.longitude = [self tanslateDegree:point.longitude];
    point.latitude = [self tanslateDegree:point.latitude];
    point.height = gaussPoint.h;
    
    return point;
}

/**
 84坐标点转80坐标点
 */
- (GaussPoint) Xian80PointFromWGS84Point:(GeodeticPoint) point {
    
    // 84坐标点先转空间直角坐标点
    RectangularPlanePoint rectPoint84 = [self rectangularPlaneCoordinateFromGeodeticCoordinates:point inEllipsoid:[Ellipsoid WGS84]];
    
    // 84空间坐标点转80空间坐标点
    MapParameter parameter;
    //杭州参数
//    parameter.deltaX = 202.634955018375;
//    parameter.deltaY = 80.880255411162;
//    parameter.deltaZ = 66.252987028366;
//    parameter.rotateX = 1.227450;
//    parameter.rotateY = 2.290472;
//    parameter.rotateZ = -2.866066;
//    parameter.k = -1.69905481167421;
    
    // 宁波参数 84转80参数
    parameter.deltaX = 219.24366659464;
    parameter.deltaY = 85.294116277735;
    parameter.deltaZ = 72.290279203019;
    parameter.rotateX = 1.150130;
    parameter.rotateY = 2.901362;
    parameter.rotateZ = -3.161284;
    parameter.k = -1.427491005181;
    
    RectangularPlanePoint rectoPoint80 = [self rectangularPlaneTanslate:rectPoint84 withParameter:parameter];
    
    // 80空间坐标点转80经纬度坐标
    GeodeticPoint geoPoint80 = [self geodeticCoordinatesFromRectangularPlaneCoordinate:rectoPoint80 inEllipsoid:[Ellipsoid Xian80]];
    
    // 80经纬度坐标点转80平面坐标
    GaussPoint toPoint = [self gaussProjection:geoPoint80 inEllipsoid:[Ellipsoid Xian80] withZoneWidth:3];
    
    return toPoint;
}

@end


@implementation Ellipsoid
@synthesize a = _a;
@synthesize b = _b;
@synthesize f = _f;
@synthesize e2 = _e2;
@synthesize ee = _ee;

- (double) ee {
    _ee = self.e2/(1 - self.e2);
    return _ee;
}

+ (Ellipsoid*) WGS84 {
    Ellipsoid* wgs84 = [[Ellipsoid alloc] init];
    wgs84.a = 6378137;
    wgs84.b = 6356752.314;
    wgs84.f = 1/298.257223563;
    wgs84.e2 = 0.006694379989;
    return wgs84;
}

+ (Ellipsoid*) Xian80 {
    Ellipsoid* xian80 = [[Ellipsoid alloc] init];
    xian80.a = 6378140;
    xian80.b = 6356755.2882;
    xian80.f = 1/298.257;
    xian80.e2 = 0.00669438499959;
    return xian80;
}

+ (Ellipsoid*) BeiJin54 {
    Ellipsoid* beijin54 = [[Ellipsoid alloc] init];
    beijin54.a = 6378245;
    beijin54.b = 6356863.019;
    beijin54.f = 1/298.3;
    beijin54.e2 = 0.006693421623;
    return beijin54;
}
@end