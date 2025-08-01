#include <cmath>

#ifndef VECTOR3D_H
#define VECTOR3D_H

struct Vector3D {
    float x, y, z;

    Vector3D(float x_=0, float y_=0, float z_=0) : x(x_), y(y_), z(z_) {}

    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }

    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }

    Vector3D operator*(float scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    float norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    float dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector3D& operator+=(const Vector3D& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    Vector3D& operator-=(const Vector3D& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    Vector3D normalized() const {
        double n = this->norm();
        if (n < 1e-12) return Vector3D{0.0, 0.0, 0.0};
        return (*this) * (1.0 / n);
    }

    Vector3D apply_minimum_image(const Vector3D& dr, const Vector3D& box_size) {
        Vector3D r = dr;
        if (box_size.x > 0) r.x -= box_size.x * std::round(r.x / box_size.x);
        if (box_size.y > 0) r.y -= box_size.y * std::round(r.y / box_size.y);
        if (box_size.z > 0) r.z -= box_size.z * std::round(r.z / box_size.z);
        return r;
    }

};

inline Vector3D apply_minimum_image(const Vector3D& dr, const Vector3D& box_size) {
    Vector3D r = dr;
    if (box_size.x > 0) r.x -= box_size.x * std::round(r.x / box_size.x);
    if (box_size.y > 0) r.y -= box_size.y * std::round(r.y / box_size.y);
    if (box_size.z > 0) r.z -= box_size.z * std::round(r.z / box_size.z);
    return r;
}

inline Vector3D wrap_position(const Vector3D& pos, const Vector3D& box_size) {
    return Vector3D(
        pos.x - std::floor(pos.x / box_size.x) * box_size.x,
        pos.y - std::floor(pos.y / box_size.y) * box_size.y,
        pos.z - std::floor(pos.z / box_size.z) * box_size.z
    );
}

#endif

