// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../aterms/cache.h"

BOOST_AUTO_TEST_SUITE(aterm_cache)

using everybeam::aterms::Cache;

BOOST_AUTO_TEST_CASE(construction) {
  Cache cache(100);

  BOOST_CHECK_EQUAL(cache.ATermSize(), 100);
  BOOST_CHECK_EQUAL(cache.Find(0), Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(150.0e6), Cache::kNotFound);
}

BOOST_AUTO_TEST_CASE(store) {
  const size_t n = 100;
  std::vector<std::complex<float>> dataA(n, 1.0), dataB(n, -3.0),
      scratch(n, 0.0);
  dataA[37] = 5.0;
  dataB[n - 1] = -5.0;

  Cache cache(n);
  cache.Store(100e6, dataA.data());
  size_t index = cache.Find(100e6);
  BOOST_CHECK_NE(index, Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(0.0), Cache::kNotFound);

  cache.Get(index, scratch.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(dataA.begin(), dataA.end(), scratch.begin(),
                                scratch.end());

  cache.Store(200e6, dataB.data());
  BOOST_CHECK_NE(cache.Find(200e6), Cache::kNotFound);
  BOOST_CHECK_NE(cache.Find(100e6), Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(0.0), Cache::kNotFound);

  cache.Get(cache.Find(200e6), scratch.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(dataB.begin(), dataB.end(), scratch.begin(),
                                scratch.end());
  cache.Get(cache.Find(100e6), scratch.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(dataA.begin(), dataA.end(), scratch.begin(),
                                scratch.end());

  cache.Store(150e6, dataA.data());
  BOOST_CHECK_NE(cache.Find(200e6), Cache::kNotFound);
  BOOST_CHECK_NE(cache.Find(150e6), Cache::kNotFound);
  BOOST_CHECK_NE(cache.Find(100e6), Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(0.0), Cache::kNotFound);
}

BOOST_AUTO_TEST_CASE(overwrite) {
  const size_t n = 100;
  std::vector<std::complex<float>> dataA(n, 1.0), dataB(n, -3.0), dataC(n, 9.0),
      scratch(n, 0.0);
  dataA[37] = 5.0;
  dataB[n - 1] = -5.0;
  dataC[n / 2] = 9.9;

  Cache cache(n);
  cache.Store(100e6, dataA.data());
  cache.Store(200e6, dataB.data());
  cache.Store(100e6, dataC.data());
  BOOST_CHECK_NE(cache.Find(200e6), Cache::kNotFound);
  BOOST_CHECK_NE(cache.Find(100e6), Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(0.0), Cache::kNotFound);

  cache.Get(cache.Find(100e6), scratch.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(dataC.begin(), dataC.end(), scratch.begin(),
                                scratch.end());
  cache.Get(cache.Find(200e6), scratch.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(dataB.begin(), dataB.end(), scratch.begin(),
                                scratch.end());
}

BOOST_AUTO_TEST_CASE(reset) {
  const size_t n = 100;
  Cache cache(n);
  cache.Reset();

  BOOST_CHECK_EQUAL(cache.ATermSize(), 100);
  BOOST_CHECK_EQUAL(cache.Find(0.0), Cache::kNotFound);

  std::vector<std::complex<float>> dataA(n, 1.0), dataB(n, -3.0), dataC(n, 9.0),
      scratch(n, 0.0);
  dataA[37] = 5.0;
  dataB[n - 1] = -5.0;
  dataC[n / 2] = 9.9;

  cache.Store(100e6, dataA.data());
  cache.Store(200e6, dataB.data());

  cache.Reset();
  BOOST_CHECK_EQUAL(cache.Find(200e6), Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(100e6), Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(0.0), Cache::kNotFound);

  cache.Store(150e6, dataB.data());
  cache.Store(200e6, dataC.data());

  BOOST_CHECK_EQUAL(cache.Find(100e6), Cache::kNotFound);
  BOOST_CHECK_EQUAL(cache.Find(0.0), Cache::kNotFound);
  BOOST_CHECK_NE(cache.Find(200e6), Cache::kNotFound);
  BOOST_CHECK_NE(cache.Find(150e6), Cache::kNotFound);

  cache.Get(cache.Find(150e6), scratch.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(dataB.begin(), dataB.end(), scratch.begin(),
                                scratch.end());
  cache.Get(cache.Find(200e6), scratch.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(dataC.begin(), dataC.end(), scratch.begin(),
                                scratch.end());
}

BOOST_AUTO_TEST_SUITE_END()
