/**
 * @file label.h
 * @brief
 * @version 0.1
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

template <class weight_t>
struct label {
    int vertex, prev;
    weight_t dist;

    label(int _v, weight_t _w, int _prev) : vertex(_v), prev(_prev), dist(_w) {}

    bool operator< (const label &u) { return vertex < u.vertex; }
};
