#include <ITRFDirection.h>
#include <ITRFConverter.h>
#include <ElementResponse.h>
#include <Station.h>
#include <LofarMetaDataUtil.h>

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/tables/Tables/TableKeyword.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>

using matrix22c_t = LOFAR::StationResponse::matrix22c_t;
using vector3r_t = LOFAR::StationResponse::vector3r_t;

void GetPhaseCentreInfo(
    casacore::MeasurementSet& ms,
    size_t fieldId,
    double& ra,
    double& dec);

void GetITRFDirections(
    vector3r_t* itrfDirections,
    size_t subgrid_size,
    float image_size,
    double time,
    double phaseCentreRA,
    double phaseCentreDec);

void StoreATermsReal(
    const std::string& filename,
    const std::complex<float>* buffer,
    size_t nStations,
    size_t width,
    size_t height);

void StoreBeam(
    const std::string& filename,
    const std::complex<float>* buffer,
    size_t nStations,
    size_t width,
    size_t height);

void setITRFVector(
    const casacore::MDirection& itrfDir,
    vector3r_t& itrf);

inline matrix22c_t operator*(
    const matrix22c_t &arg0,
    const matrix22c_t &arg1)
{
    matrix22c_t result;
    result[0][0] = arg0[0][0] * arg1[0][0] + arg0[0][1] * arg1[1][0];
    result[0][1] = arg0[0][0] * arg1[0][1] + arg0[0][1] * arg1[1][1];
    result[1][0] = arg0[1][0] * arg1[0][0] + arg0[1][1] * arg1[1][0];
    result[1][1] = arg0[1][0] * arg1[0][1] + arg0[1][1] * arg1[1][1];
    return result;
}

template<typename T>
void XYToLM(
    size_t x,
    size_t y,
    T pixelSizeX,
    T pixelSizeY,
    size_t width,
    size_t height,
    T &l,
    T &m)
{
	T midX = (T) width / 2.0, midY = (T) height / 2.0;
	l = (midX - (T) x) * pixelSizeX;
	m = ((T) y - midY) * pixelSizeY;
}