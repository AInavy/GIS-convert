#include <iostream>
#include "proj_api.h"
#include <iomanip>
using namespace std;

void ecef_to_lla(double x, double y, double z, double &latitude, double &longitude, double &alt);
void lla_to_ecef(double latitude, double longitude, double alt, double &x, double &y, double &z);
void WGS842GKP(double Latitude, double Longitude, double Hell, double &x, double &y, double &z);
void GKP2WGS84(double x, double y, double z, int zoneID, double &Longitude, double &Latitude, double &Hell);

#define  PI  3.141592653
double WGS84_A = 6378137.0;
double WGS84_B = 6356752.314245;
double WGS84_E = 0.0818191908;

void lla_to_ecef(double latitude, double longitude, double alt, double &x, double &y, double &z)
{
    double clat = cos(latitude*PI / 180);
    double slat = sin(latitude*PI / 180);
    double clon = cos(longitude*PI / 180);
    double slon = sin(longitude*PI / 180);

    double a2 = WGS84_A * WGS84_A;
    double b2 = WGS84_B * WGS84_B;

    double L = 1.0 / sqrt(a2 * clat*clat + b2 * slat*slat);
    x = (a2 * L + alt) * clat * clon;
    y = (a2 * L + alt) * clat * slon;
    z = (b2 * L + alt) * slat;

}

void ecef_to_lla(double x, double y, double z, double &latitude, double &longitude, double &alt)
{
    double b = sqrt(WGS84_A*WGS84_A*(1 - WGS84_E * WGS84_E));
    double ep = sqrt((WGS84_A*WGS84_A - b * b) / (b*b));
    double p = hypot(x, y);
    double th = atan2(WGS84_A*z, b*p);
    double lon = atan2(y, x);
    double lat = atan2((z + ep * ep*b* pow(sin(th), 3)), (p - WGS84_E * WGS84_E*WGS84_A*pow(cos(th), 3)));
    double N = WGS84_A / sqrt(1 - WGS84_E * WGS84_E*sin(lat)*sin(lat));
    double alt_ = p / cos(lat) - N;

    latitude = lat / PI * 180.0;
    longitude = lon / PI * 180.0;
    alt = alt_;

}

void WGS842GKP(double Latitude, double Longitude, double Hell, double &x, double &y, double &z)
{
    projPJ pj_merc, pj_latlong;

    char pjValue[256] = { 0 };

    int zoneID = int((Longitude + 1.5) / 3);

    sprintf(pjValue, "+proj=tmerc +lon_0=%d +k=1 +x_0=500000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", (zoneID * 3));

    if (!(pj_merc = pj_init_plus(pjValue)))
        exit(1);
    if (!(pj_latlong = pj_init_plus("+proj=latlong +ellps=WGS84 +datum=WGS84")))
        exit(1);

    x = Longitude;
    y = Latitude;
    z = Hell;

    x *= DEG_TO_RAD;
    y *= DEG_TO_RAD;


    pj_transform(pj_latlong, pj_merc, 1, 1, &x, &y, &z);

    pj_free(pj_merc);
    pj_free(pj_latlong);

}

void GKP2WGS84(double x, double y, double z, int zoneID, double &Longitude, double &Latitude, double &Hell)
{
    projPJ pj_merc, pj_latlong;
    char pjValue[256] = { 0 };

    sprintf(pjValue, "+proj=tmerc +lon_0=%d +k=1 +x_0=500000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", zoneID * 3);

    if (!(pj_merc = pj_init_plus(pjValue)))
        exit(1);
    if (!(pj_latlong = pj_init_plus("+proj=latlong +ellps=WGS84 +datum=WGS84")))
        exit(1);



    pj_transform(pj_merc, pj_latlong, 1, 1, &x, &y, &z);

    x *= RAD_TO_DEG;
    y *= RAD_TO_DEG;


    Longitude = x;
    Latitude = y;
    Hell = z;

    pj_free(pj_merc);
    pj_free(pj_latlong);

}


int main()
{

    double lon = 115.50000000; //116.24110390;
    double lat = 40.07124398; //40.06720701;
    double x = 0.0, y = 0.0, z = 0.0;
    WGS842GKP(lat, lon, 0.0, x, y, z);
    cout<<"wgs84:           "<<setiosflags(ios::fixed)<<setprecision(8)<<lon<<", "<<lat<<endl;
    cout<<"converted gauss: "<<setiosflags(ios::fixed)<<setprecision(8)<<x<<", "<<y<<endl;
    double lon1 = 115.49999999; //116.24110390;
    double lat1 = 40.07124398; //40.06720701;
    WGS842GKP(lat1, lon1, 0.0, x, y, z);
    cout<<"wgs84:           "<<setiosflags(ios::fixed)<<setprecision(8)<<lon1<<", "<<lat1<<endl;
    cout<<"converted gauss: "<<setiosflags(ios::fixed)<<setprecision(8)<<x<<", "<<y<<endl;
    return 0;
}
