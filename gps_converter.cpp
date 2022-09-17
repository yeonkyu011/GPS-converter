#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>


using namespace std;
using namespace Eigen;

vector<double> geodetic2ECEF(double latitude,double longitude,double altitude);
double deg2Rad(const double degrees);
Matrix3d nRe(double latitude, double longitude);
Vector3d ECEF2ENU(double x0, double y0, double z0, double x, double y, double z, Matrix3d ECEF2Ned_Mat);



static double kSemimajorAxis = 6378137;
static double kSemiminorAxis = 6356752.3142;
static double kFirstEccentricitySquared = 6.69437999014 * 0.001;
static double kSecondEccentricitySquared = 6.73949674228 * 0.001;
static double kFlattening = 1 / 298.257223563;

double deg2Rad(const double degrees)
{
    return (degrees/180.0)*M_PI;
}

vector<double> geodetic2ECEF(double latitude, double longitude, double altitude)
{
    double lat_rad = deg2Rad(latitude);
    double lon_rad = deg2Rad(longitude);
    double xi = sqrt(1 - kFirstEccentricitySquared * sin(lat_rad) * sin(lat_rad));
    double x = (kSemimajorAxis / xi + altitude) * cos(lat_rad) * cos(lon_rad);
    double y = (kSemimajorAxis / xi + altitude) * cos(lat_rad) * sin(lon_rad);
    double z = (kSemimajorAxis / xi * (1 - kFirstEccentricitySquared) + altitude) * sin(lat_rad);
    vector<double> ecef_array = {x,y,z};
    return ecef_array;
}

Matrix3d nRe(double latitude, double longitude)
{
    double rad_lat = deg2Rad(latitude);
    double rad_lon = deg2Rad(longitude);    
    double sLat = sin(rad_lat);
    double sLon = sin(rad_lon);
    double cLat = cos(rad_lat);
    double cLon = cos(rad_lon);

    Matrix3d ret;
    ret(0, 0) = -sLon;
    ret(0, 1) = -sLat * cLon;
    ret(0, 2) = cLat*cLon;
    ret(1, 0) = cLon;
    ret(1, 1) = -sLat*sLon;
    ret(1, 2) = cLat*sLon;
    ret(2, 0) = 0;
    ret(2, 1) = cLat;
    ret(2, 2) = sLat;

    Matrix3d enu;
    enu = ret.transpose();

    return enu;
}

Vector3d ECEF2ENU(double x0, double y0, double z0, double x, double y, double z, Matrix3d ECEF2ENU_Mat)
{
    Vector3d ENU_array0;
    ENU_array0[0] = x - x0;
    ENU_array0[1] = y - y0;
    ENU_array0[2] = z - z0;
    Vector3d ENU_array;
    ENU_array = ECEF2ENU_Mat * ENU_array0;
    return ENU_array; // north,east,up
}

Vector2d xytransform(double angle, double North, double East)
{
    Matrix2d xytransform;
    xytransform(0,0) = cos(angle);
    xytransform(0,1) = -sin(angle);
    xytransform(1,0) = sin(angle);
    xytransform(1,1) = cos(angle);

    Vector2d ENarray;
    ENarray[0] = North;
    ENarray[1] = East;

    Vector2d xyarray;
    xyarray = xytransform * ENarray;

    xyarray[1] = -xyarray[1];

    return xyarray;   
}

int main()
{
    vector<double> initial(3);
    vector<double> destination(3);
    double angle;

    cout << "angle: "; cin >>angle;
    cout << "Initial latitude: " ; cin >>initial[0];
    cout << "Initial longitude: " ; cin >>initial[1];
    cout << "Initial altitude: " ; cin >>initial[2];
    
    angle = 360-angle;
    angle = deg2Rad(angle);

    cout << "Destination latitude: " ; cin >>destination[0];
    cout << "Destination longitude: " ; cin >>destination[1];
    cout << "Destination altitude: " ; cin >>destination[2];

    vector<double> ECEF_init(3);
    vector<double> ECEF_dest(3);
    ECEF_init = geodetic2ECEF(initial[0],initial[1],initial[2]);
    ECEF_dest = geodetic2ECEF(destination[0],destination[1],destination[2]);
    Matrix3d ECEF2ENU_Mat;
    ECEF2ENU_Mat = nRe(initial[0], initial[1]);
    Vector3d ENU;
    ENU = ECEF2ENU(ECEF_init[0],ECEF_init[1],ECEF_init[2],ECEF_dest[0],ECEF_dest[1],ECEF_dest[2],ECEF2ENU_Mat);
    cout << "EAST: " << ENU[0] << " NORTH: " << ENU[1] << " UP: " << ENU[2] << endl;
    Vector2d xyarray;
    xyarray = xytransform(angle, ENU[1], ENU[0]);
    cout << "X: " << xyarray[0] << " Y: " << xyarray[1] << endl;

    return 0;
}
