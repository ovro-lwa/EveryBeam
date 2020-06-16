#ifndef EVERYBEAM_BEAMFORMER_H
#define EVERYBEAM_BEAMFORMER_H

#include <complex>
#include <vector>

#include "Element.h"
#include "Types.h"

namespace everybeam {

class BeamFormer : public Antenna
{
public:

    typedef std::shared_ptr<BeamFormer> Ptr;

    BeamFormer() :
        Antenna()
    {
        m_local_phase_reference_position = transform_to_local_position(m_phase_reference_position);
    }

    BeamFormer(const CoordinateSystem &coordinate_system) :
        Antenna(coordinate_system)
    {
        m_local_phase_reference_position = transform_to_local_position(m_phase_reference_position);
    }

    BeamFormer(CoordinateSystem coordinate_system, vector3r_t phase_reference_position) :
        Antenna(coordinate_system, phase_reference_position)
    {
        m_local_phase_reference_position = transform_to_local_position(m_phase_reference_position);
    }

    void add_antenna(Antenna::Ptr antenna) {m_antennas.push_back(antenna);}

private:

    vector3r_t  m_local_phase_reference_position; // in coordinate system of Antenna

    vector3r_t transform_to_local_position(const vector3r_t &position);

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

    std::vector<std::complex<double>> compute_geometric_response(double freq, const vector3r_t &direction) const;
    std::vector<std::pair<std::complex<double>,std::complex<double>>> compute_weights(const vector3r_t &direction, double freq) const;

    std::vector<Antenna::Ptr> m_antennas;


};

} // namespace everybeam


#endif
