/**
 * @file label.h
 * @brief
 * @version 0.1
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <iostream>

template <class weight_t>
struct label {
    size_t vertex, prev;
    weight_t dist;

    label(size_t _v, weight_t _w, size_t _prev) : vertex(_v), prev(_prev), dist(_w) {}

    bool operator<(const label &u) { return vertex < u.vertex; }
};

template <class weight_t>
std::ostream &operator<<(std::ostream &os, const label<weight_t> &cur) {
    return os << "(" << cur.vertex << ", " << cur.dist << ", " << cur.prev << ")";
}
