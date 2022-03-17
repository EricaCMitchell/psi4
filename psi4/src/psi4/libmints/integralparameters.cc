/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/integralparameters.h"

using namespace psi;

CorrelationFactor::CorrelationFactor(size_t nparam) : IntegralParameters(nparam) {}

CorrelationFactor::CorrelationFactor(std::shared_ptr<Vector> coeff, std::shared_ptr<Vector> exponent)
    : IntegralParameters(coeff->dim()) {
    set_params(coeff, exponent);
}

CorrelationFactor::~CorrelationFactor() {
    delete[] coeff_;
    delete[] exponent_;
}

void CorrelationFactor::set_params(std::shared_ptr<Vector> coeff, std::shared_ptr<Vector> exponent) {
    int nparam = coeff->dim();
    if (nparam) {
        coeff_ = new double[nparam];
        exponent_ = new double[nparam];
        for (int i = 0; i < nparam; ++i) {
            coeff_[i] = coeff->get(0, i);
            exponent_[i] = exponent->get(0, i);
        }
    }
}

FittedSlaterCorrelationFactor::FittedSlaterCorrelationFactor(double exponent) : CorrelationFactor(6) {
    // Perform the fit.
    auto exps = std::make_shared<Vector>(6);
    auto coeffs = std::make_shared<Vector>(6);

    slater_exponent_ = exponent;

    // The fitting coefficients
    coeffs->set(0, 0, -0.31442480597241274);
    coeffs->set(0, 1, -0.30369575353387201);
    coeffs->set(0, 2, -0.16806968430232927);
    coeffs->set(0, 3, -0.098115812152857612);
    coeffs->set(0, 4, -0.060246640234342785);
    coeffs->set(0, 5, -0.037263541968504843);

    // and the exponents
    exps->set(0, 0, 0.22085085450735284);
    exps->set(0, 1, 1.0040191632019282);
    exps->set(0, 2, 3.6212173098378728);
    exps->set(0, 3, 12.162483236221904);
    exps->set(0, 4, 45.855332448029337);
    exps->set(0, 5, 254.23460688554644);

    // They just need to be scaled
    double expsq = exponent * exponent;
    exps->scale(expsq);
    set_params(coeffs, exps);
}
