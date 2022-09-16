// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../elementresponsefixeddirection.h"

#include <boost/test/unit_test.hpp>

#include "../common/mathutils.h"

using everybeam::ElementResponse;
using everybeam::ElementResponseFixedDirection;
using everybeam::ElementResponseModel;

namespace {
const ElementResponseModel kModel = ElementResponseModel::kSkaMidAnalytical;
const int kElementId = 42;
const double kFrequency = 42.0e6;
const everybeam::vector3r_t kDirection{0.42, 1.42, 2.42};
const everybeam::vector2r_t kThetaPhi = everybeam::cart2thetaphi(kDirection);
const double kTheta = kThetaPhi[0];
const double kPhi = kThetaPhi[1];
const aocommon::MC2x2 kResponse{{1, 2}, {3, 4}, {5, 6}, {7, 8}};

// Give these values to ElementResponseFixedDirection and check in
// MockElementResponse that it receives the fixated values instead.
const double kWrongTheta = 0.0;
const double kWrongPhi = 0.0;

class MockElementResponse : public everybeam::ElementResponse {
 public:
  MockElementResponse() {}

  ElementResponseModel GetModel() const override { return kModel; }

  aocommon::MC2x2 Response(double frequency, double theta,
                           double phi) const override {
    BOOST_TEST(frequency == kFrequency);
    BOOST_TEST(theta == kTheta);
    BOOST_TEST(phi == kPhi);
    ++response_without_id_count_;
    return kResponse;
  }

  aocommon::MC2x2 Response(int element_id, double frequency, double theta,
                           double phi) const override {
    BOOST_TEST(element_id == kElementId);
    BOOST_TEST(frequency == kFrequency);
    BOOST_TEST(theta == kTheta);
    BOOST_TEST(phi == kPhi);
    ++response_with_id_count_;
    return kResponse;
  }

  mutable int response_without_id_count_ = 0;
  mutable int response_with_id_count_ = 0;
};

void CheckResponses(const ElementResponse& element_response,
                    const MockElementResponse& mock) {
  const aocommon::MC2x2 response_without_id =
      element_response.Response(kFrequency, kWrongTheta, kWrongPhi);
  BOOST_TEST(mock.response_without_id_count_ == 1);
  BOOST_TEST(mock.response_with_id_count_ == 0);
  BOOST_CHECK_EQUAL_COLLECTIONS(response_without_id.Data(),
                                response_without_id.Data() + 4,
                                kResponse.Data(), kResponse.Data() + 4);

  const aocommon::MC2x2 response_with_id =
      element_response.Response(kElementId, kFrequency, kWrongTheta, kWrongPhi);
  BOOST_TEST(mock.response_without_id_count_ == 1);
  BOOST_TEST(mock.response_with_id_count_ == 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(response_with_id.Data(),
                                response_with_id.Data() + 4, kResponse.Data(),
                                kResponse.Data() + 4);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(element_response_fixed_direction)

BOOST_AUTO_TEST_CASE(get_model) {
  auto mock = std::make_shared<MockElementResponse>();
  ElementResponseFixedDirection direct(mock, 0.0, 0.0);
  BOOST_TEST(direct.GetModel() == kModel);

  std::shared_ptr<ElementResponse> fixated = mock->FixateDirection(kDirection);
  BOOST_TEST(dynamic_cast<ElementResponseFixedDirection*>(fixated.get()));
  BOOST_TEST(fixated->GetModel() == kModel);
}

BOOST_AUTO_TEST_CASE(response_direct) {
  // Test using a directly-created ElementResponseFixedDirection.
  auto mock = std::make_shared<MockElementResponse>();
  ElementResponseFixedDirection direct(mock, kTheta, kPhi);
  CheckResponses(direct, *mock);
}

BOOST_AUTO_TEST_CASE(response_fixated) {
  // Test using an object from ElementResponse::FixateDirection.
  auto mock = std::make_shared<MockElementResponse>();
  std::shared_ptr<ElementResponse> fixated = mock->FixateDirection(kDirection);
  CheckResponses(*fixated, *mock);
}

BOOST_AUTO_TEST_CASE(fixate_direction) {
  // Test the FixateDirection of ElementResponseFixedDirection itself.
  auto mock = std::make_shared<MockElementResponse>();
  auto fixated = std::make_shared<ElementResponseFixedDirection>(
      mock, kWrongTheta, kWrongPhi);
  std::shared_ptr<ElementResponse> fixated_twice =
      fixated->FixateDirection(kDirection);
  BOOST_CHECK(
      dynamic_cast<ElementResponseFixedDirection*>(fixated_twice.get()));

  BOOST_TEST(fixated_twice->GetModel() == kModel);
  CheckResponses(*fixated_twice, *mock);
}

BOOST_AUTO_TEST_SUITE_END()
