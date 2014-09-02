//
//  main.m
//  CoordinateUtils
//
//  Created by damingdan on 14/8/27.
//  Copyright (c) 2014年 Kingoit. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "CoordinateUtils.h"
/**
 2014-08-29 17:07:16.933 MobileMap[771:60b] update location (gauss point) x:3349948.338572, y:40509213.081722
 2014-08-29 17:07:16.934 MobileMap[771:60b] update location x:120.095745, y:30.269101
 
 
 2014-08-29 17:12:28.943 MobileMap[777:60b] update location (gauss point) x:3349943.346089, y:40509215.376527
 2014-08-29 17:12:28.944 MobileMap[777:60b] update location x:120.095769, y:30.269055
 2014-08-29 17:12:28.947 MobileMap[777:60b] update location x:509215.380868, y:3349944.907049
 
 
 update location (gauss point) x:30.269946, y:120.095285, h:8.587585
 */

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        CoordinateUtils* coordinateUtils = [[CoordinateUtils alloc] init];
        GeodeticPoint geoPoint;
        geoPoint.latitude = 29.48413400 ;
        geoPoint.longitude = 121.3621240 ;
        geoPoint.height = 15.029;
        GaussPoint point = [coordinateUtils gaussProjection:geoPoint inEllipsoid:[Ellipsoid WGS84] withZoneWidth:3];
        NSLog(@"%f, %f", point.x, point.y);
        
//        point.x = 3300304.414;
//        point.y = 655131.668;
//        point.h = 1.996;
        geoPoint = [coordinateUtils gaussProjectionReverse:point inEllipsoid:[Ellipsoid WGS84] withZoneWidth:3];
        NSLog(@"%f, %f", geoPoint.latitude, geoPoint.longitude);
        
//        RectangularPlanePoint rectPoint = [coordinateUtils rectangularPlaneCoordinateFromGeodeticCoordinates:geoPoint inEllipsoid:[Ellipsoid WGS84]];
//        NSLog(@"%f, %f, %f", rectPoint.x, rectPoint.y, rectPoint.z);
//        
//        geoPoint = [coordinateUtils geodeticCoordinatesFromRectangularPlaneCoordinate:rectPoint inEllipsoid:[Ellipsoid WGS84]];
//        NSLog(@"%f, %f, %f", geoPoint.latitude, geoPoint.longitude, geoPoint.height);
        
    }
    
    return 0;
}
