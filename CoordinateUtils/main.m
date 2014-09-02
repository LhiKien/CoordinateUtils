//
//  main.m
//  CoordinateUtils
//
//  Created by damingdan on 14/8/27.
//  Copyright (c) 2014年 Kingoit. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "CoordinateUtils.h"

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        CoordinateUtils* coordinateUtils = [[CoordinateUtils alloc] init];
        GeodeticPoint geoPoint;
        geoPoint.latitude = 29.48413400 ;
        geoPoint.longitude = 121.3621240 ;
        geoPoint.height = 15.029;
//        GaussPoint point = [coordinateUtils gaussProjection:geoPoint inEllipsoid:[Ellipsoid WGS84] withZoneWidth:3];
//        NSLog(@"%f, %f", point.x, point.y);
//        
//        geoPoint = [coordinateUtils gaussProjectionReverse:point inEllipsoid:[Ellipsoid WGS84] withZoneWidth:3];
//        NSLog(@"%f, %f", geoPoint.latitude, geoPoint.longitude);
//        
//        RectangularPlanePoint rectPoint = [coordinateUtils rectangularPlaneCoordinateFromGeodeticCoordinates:geoPoint inEllipsoid:[Ellipsoid WGS84]];
//        NSLog(@"%f, %f, %f", rectPoint.x, rectPoint.y, rectPoint.z);
//        
//        geoPoint = [coordinateUtils geodeticCoordinatesFromRectangularPlaneCoordinate:rectPoint inEllipsoid:[Ellipsoid WGS84]];
//        NSLog(@"%f, %f, %f", geoPoint.latitude, geoPoint.longitude, geoPoint.height);
        
        // 宁波参数 84转80参数
        MapParameter parameter;
        parameter.deltaX = 219.24366659464;
        parameter.deltaY = 85.294116277735;
        parameter.deltaZ = 72.290279203019;
        parameter.rotateX = 1.150130;
        parameter.rotateY = 2.901362;
        parameter.rotateZ = -3.161284;
        parameter.k = -1.427491005181;
        
        RectangularPlanePoint rectPoint84;
        rectPoint84.x = -2902680.477;
        rectPoint84.y = 4717151.300;
        rectPoint84.z = 3152266.518;
        
        
        GaussPoint point = [coordinateUtils Xian80PointFromWGS84Point:geoPoint];
        
        NSLog(@"%f, %f, %f", point.x, point.y, point.h);
    }
    
    return 0;
}
