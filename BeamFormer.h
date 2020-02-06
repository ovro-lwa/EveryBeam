#ifndef BEAMFORMER_H
#define BEAMFORMER_H

#include <complex>
#include <vector>

#include "Element.h"
#include "Types.h"

namespace LOFAR {
namespace StationResponse {

class BeamFormer : public Antenna
{
public:

    typedef std::shared_ptr<BeamFormer> Ptr;

    void add_antenna(Antenna::Ptr antenna) {m_antennas.push_back(antenna);}

private:

    virtual matrix22c_t local_response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const override;

    virtual diag22c_t local_arrayFactor(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const override
    {
        return {1.0, 1.0};
    }

    vector3r_t compute_local_pointing(double time) const;
    std::vector<std::complex<double>> compute_geometric_response(double freq, const vector3r_t &direction) const;
    std::vector<std::pair<std::complex<double>,std::complex<double>>> compute_weights(double time, double freq) const;

    std::vector<Antenna::Ptr> m_antennas;

    vector2r_t m_pointing;  // ra, dec

};

} // namespace StationResponse
} // namespace LOFAR


#endif
