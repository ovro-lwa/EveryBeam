#ifndef ELEMENT_H
#define ELEMENT_H

#include <complex>
#include <memory>

#include "Antenna.h"
#include "ElementResponse.h"
#include "Types.h"

namespace LOFAR {
namespace StationResponse {

class Element : public Antenna
{
public:

    typedef std::shared_ptr<Element> Ptr;

    Element(const CoordinateSystem &coordinate_system, ElementResponse::Ptr element_response, int id) :
        Antenna(coordinate_system),
        m_id(id),
        m_element_response(element_response)
    {}

    matrix22c_t local_response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        size_t id,
        const Options &options) const;

private:

    virtual matrix22c_t local_response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const final override;

    int m_id;
    ElementResponse::Ptr m_element_response;
};

} // namespace StationResponse
} // namespace LOFAR

#endif