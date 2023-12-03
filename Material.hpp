//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType
{
    DIFFUSE,
    MICROFACET,
    REFRACT,
    REFLECT
};

class Material{
private:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // 计算折射方向
    Vector3f refract2(const Vector3f &wi,const Vector3f &N, const float &ior)
    {
        // 计算入射角余弦
        float cosi = dotProduct(wi, N);
        float eta = ior;
        // 判断是否从内向外折射
        bool out = cosi > 0;
        Vector3f n = N;
        // 如果是，反转法线方向和折射率比例
        if (out)
        {
            n = -n;
            eta = 1 / eta;
            cosi = -cosi;
        }
        // 计算判别式
        float k = 1 - eta * eta * (1 - cosi * cosi);
        // 如果判别式小于零，发生全反射，返回零向量
        if (k < 0)
            return Vector3f(0, 0, 0);
        // 否则，根据公式计算折射方向并返回
        return wi * eta + n * (eta * cosi - sqrtf(k));
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) // 光线在外部，
        { 
            cosi = -cosi; 
        }  
        else// 光线在内部
        { 
            std::swap(etai, etat); n= -N;
        }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi); // 小于0 -- 全反射 ； 大于等于0 -- 有折射
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

    float D_F(float NdotH, float roughness)
    {
        float a = roughness * roughness;
        float a2 = a * a;
        float m = NdotH * NdotH * (a2 - 1) + 1;
        return a2 / std::max((M_PI * m * m), 0.000001f);
        // return a2 / (M_PI * m * m);
    }

    float SmithG_GGX(float NdotV, float roughness)
    {
        float r = 0.5 + roughness / 2.0f;
        float m = r * r + (1 - r * r) * NdotV * NdotV;
        return 2.0f * NdotV / (NdotV + std::sqrt(m));
    }

    // UE4 Way
    float GeometrySchlickGGX(float NdotV , float roughness)
    {
        float r = roughness + 1;
        float k = r * r / 8;
        float m = NdotV / NdotV * (1.0f - k) + k;

        return NdotV / m; 
    }

    float G_F(Vector3f N, Vector3f V, Vector3f L, float roughness)
    {
        float NdotV = std::max(dotProduct(N, V), 0.0f);
        float NdotL = std::max(dotProduct(N, L), 0.0f);
        // float ggx1 = SmithG_GGX(NdotL, roughness);
        // float ggx2 = SmithG_GGX(NdotV, roughness);
        float ggx1 = GeometrySchlickGGX(NdotL, roughness);
        float ggx2 = GeometrySchlickGGX(NdotV, roughness);
        return ggx1 * ggx2;
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    float roughness;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) 
{
    return Vector3f();
}

Vector3f Material::sample(const Vector3f &wi, const Vector3f &N)
{
    switch(m_type)
    {
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);
            break;
        }
        case MICROFACET:
        {
            // iuniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r * std::cos(phi), r * std::sin(phi), z);
            return toWorld(localRay, N);
            break;
        }
        case REFLECT:
        {
            // importance sample on the hemisphere
            return (2 * normalize(N) * dotProduct(-wi, N) + wi).normalized(); // 完全镜面反射
            break;
        }
        case REFRACT:
        {
            // importance sample on the hemisphere
            // printf("%f\n",ior);
            // return wi - N * 2 * dotProduct(N,wi);
            return refract(wi, N, ior); // 求折射方向
            break;
        }
    }
}

// 均匀分布采样计算
float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N)
{
    switch(m_type)
    {
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case MICROFACET:
        {
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case REFLECT: // 完全反射
        {
            if (dotProduct(wo, N) > 0.0f)
            {
                // 计算反射的一致性，重要性采样
                // if ((-wi + wo - 2 * N * dotProduct(N, wo)).norm() < 0.01f) // 出射为镜面反射方向
                // {
                //     return 1.0f;
                // }
                return 1.0f;
            }
            else
                return 0.0f;
            break;
        }
        case REFRACT: // 完全折射
        {
            if (dotProduct(wo, N) > 0.0f)
                return 1.0f;
            else
                return 0.0f;
            break;
        }

    }
}

// 求物体材质表面的反射率
Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N)
{
    switch(m_type)
    {
        case DIFFUSE:
        {
            // calculate the contribution of diffuse model
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) 
            {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MICROFACET:
        {
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f)
            {
                // Cook-Torrance BRDF
                // 菲涅尔项，反射比例
                float F;
                // std::cout << ior << std::endl;

                // 计算反射部分占比，F
                fresnel(wi, N, ior, F);

                // 折射比例（即漫反射部分）
                float k_d = 1 - F;

                // 反射比例，不是颜色
                float k_s = F;

                Vector3f f_lambert, f_cook_torrance;

                // 折射，即漫反射项
                f_lambert = Kd / M_PI;

                // 反射项
                // 半程向量
                // float roughness = 0.01f;
                Vector3f H = normalize(-wi + wo);
                float NdotH = std::max(dotProduct(N, H), 0.0f);
                float NdotV = std::max(dotProduct(N, wo), 0.0f);
                float NdotL = std::max(dotProduct(N, -wi), 0.0f);

                // 法线分布函数
                float D = D_F(NdotH, roughness);

                // 阴影遮挡函数
                float G = G_F(N, wo, -wi, roughness);

                float m = 4 * std::max(NdotL * NdotV, 0.00001f);

                f_cook_torrance = G * F * D / m;
                // printf("%f", Ks );
                // std::cout << "D" << D << std::endl;
                // std::cout << "G" << G << std::endl;
                // std::cout << "m" << m << std::endl;
                // Vector3f result = k_d * f_lambert + Ks * f_cook_torrance;
                // std::cout << "result" << result << std::endl;
                return k_d * f_lambert + Ks * f_cook_torrance; //  k_d是比例，Ks是镜面反射颜色，金属为Vector3f，非金属为float
                // return Ks * f_cook_torrance;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case REFLECT: // 完全反射
        {
            if ((-wi + wo - 2 * N * dotProduct(N, wo)).norm() < 0.1f) // 出入射光线完全镜面反射才计算，否则贡献为0 cosalpha1 - cosalpha2 < 0.1f
            {
                return Vector3f(1.0f, 1.0f, 1.0f);
            }
            return Vector3f(0.0f,0.0f,0.0f);
            break;
        }
        case REFRACT: // 完全折射
        {
            if (1.0f - dotProduct(wo,refract(wi,N,ior).normalized()) <= 0.002f) // 出入射光线完全镜面反射才计算，否则贡献为0 cosalpha1 - cosalpha2 < 0.1f
            {
                // printf("Light \n");
                return Vector3f(1.0f, 1.0f, 1.0f);
            }
            return Vector3f(0.0f, 0.0f, 0.0f);
            break;
        }
        }
}

#endif //RAYTRACING_MATERIAL_H
