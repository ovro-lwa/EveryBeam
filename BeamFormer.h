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

    /**
     * @brief Construct a new BeamFormer object
     * 
     */    
    BeamFormer() :
        Antenna()
    {
        m_local_phase_reference_position = transform_to_local_position(m_phase_reference_position);
    }

    /**
     * @brief Construct a new BeamFormer object given a coordinate system.
     * 
     * @param coordinate_system 
     */    
    BeamFormer(const CoordinateSystem &coordinate_system) :
        Antenna(coordinate_system)
    {
        m_local_phase_reference_position = transform_to_local_position(m_phase_reference_position);
    }

    /**
     * @brief Construct a new BeamFormer object given a coordinate system and a phase reference position
     * 
     * @param coordinate_system 
     * @param phase_reference_position 
     */
    BeamFormer(CoordinateSystem coordinate_system, vector3r_t phase_reference_position) :
        Antenna(coordinate_system, phase_reference_position)
    {
        m_local_phase_reference_position = transform_to_local_position(m_phase_reference_position);
    }

    /**
     * @brief Add an antenna to the m_antenna array.
     * 
     * @param antenna 
     */
    void add_antenna(Antenna::Ptr antenna) {m_antennas.push_back(antenna);}

private:

    vector3r_t  m_local_phase_reference_position; // in coordinate system of Antenna

    // Transform position vector into a local position vector
    vector3r_t transform_to_local_position(const vector3r_t &position);

    // Compute the BeamFormer response in certain direction of arrival (ITRF, m)
    // and return (Jones) matrix of response
    virtual matrix22c_t local_response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const override;

    // Compute the local arrayFactor, with arrayFactor a vectorial "representation"
    // of Jones matrix
    virtual diag22c_t local_arrayFactor(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const override
    {
        return {1.0, 1.0};
    }
    
    // Compute the geometric response for all the antennas in the BeamFormer based on 
    // the probing frequency and a specified direction (either pointing dir or dir of interest).
    std::vector<std::complex<double>> compute_geometric_response(const double freq, const vector3r_t &direction) const;
    
    // Compute the weights based on the pointing direction of the beam and the beam reference frequence.
    std::vector<std::pair<std::complex<double>,std::complex<double>>> compute_weights(const vector3r_t &direction, double freq) const;

    // List of antennas in BeamFormer
    // TODO: Maybe refactor to _m_antennas to indicate m_antennas is a private attribute
    std::vector<Antenna::Ptr> m_antennas;
};
} // namespace everybeam
#endif
