#ifndef MEDYAN_Structure_SurfaceMesh_MTriangle_Hpp
#define MEDYAN_Structure_SurfaceMesh_MTriangle_Hpp

/******************************************************************************

Storing some mechanical properties of the triangle patches

******************************************************************************/

class MTriangle {

private:

    double _kExVol; ///< Local excluded volume constant, which describes
                    ///< excluded volume interactions between triangles and cylinder tips

public:
    MTriangle(short membraneType);

    void setExVolConst(double kExVol) { _kExVol = kExVol; }
    double getExVolConst() const { return _kExVol; }

};


#endif
