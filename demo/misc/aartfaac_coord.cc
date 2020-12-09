// Piece of not yet active code. When activated, it should be plugged-in
// into cpp/msreadutils.cc

// AntennaField::CoordinateSystem readCoordinateSystemAartfaac(
//     const Table &table, unsigned int id)
// {
//     ROArrayQuantColumn<Double> c_position(table, "POSITION", "m");
//
//     // Read antenna field center (ITRF).
//     Vector<Quantity> aips_position = c_position(id);
//     assert(aips_position.size() == 3);
//
//     vector3r_t position = {{aips_position(0).getValue(),
//         aips_position(1).getValue(), aips_position(2).getValue()}};
//
//     TableRecord keywordset = table.keywordSet();
//     Matrix<double> aips_axes;
//     keywordset.get("AARTFAAC_COORDINATE_AXES", aips_axes);
//     assert(aips_axes.shape().isEqual(IPosition(2, 3, 3)));
//
//     vector3r_t p = {{aips_axes(0, 0), aips_axes(1, 0), aips_axes(2, 0)}};
//     vector3r_t q = {{aips_axes(0, 1), aips_axes(1, 1), aips_axes(2, 1)}};
//     vector3r_t r = {{aips_axes(0, 2), aips_axes(1, 2), aips_axes(2, 2)}};
//
//     AntennaField::CoordinateSystem system = {position, {p, q, r}};
//
//     return system;
// }