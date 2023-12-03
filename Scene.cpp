//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include <random>

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const // 在光源位置随机采样
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) 
    {
        if (objects[k]->hasEmit())
        {
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum) // 
            {
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}
// Implementation of Path Tracing
// Vector3f Scene::castRay(const Ray &ray, int depth) const
// {
//     Intersection inter = intersect(ray);

//     if (inter.happened) // 如果光线与场景有交点
//     {
//         // 如果打到光源
//         if (inter.m->hasEmission())
//         {
//             // 如果射线第一次打到光源，直接返回光源光
//             if (depth == 0)
//             {
//                 return inter.m->getEmission();
//             }
//             // 若非射线经弹射打到光源，则在打到物体是判断，这里不做处理，返回0
//             else
//                 return Vector3f(0, 0, 0);
//         }

//         // 如果打到物体
//         Vector3f L_dir(0, 0, 0);
//         Vector3f L_indir(0, 0, 0);

//         // 均匀采样光源物体，取光源上一点
//         Intersection lightInter;
//         float pdf_light = 0.0f;
//         sampleLight(lightInter, pdf_light);

//         // 物体表面法线
//         auto &N = inter.normal;
//         // 灯光表面法线
//         auto &NN = lightInter.normal;
//         // 物体点坐标
//         auto &objPos = inter.coords;
//         // 光源点坐标
//         auto &lightPos = lightInter.coords;

//         auto diff = lightPos - objPos;
//         auto lightDir = diff.normalized();
//         float lightDistance = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;

//         // 从物体打向光源的感光线(感光线为光线传播的逆方向)
//         Ray light(objPos, lightDir);
//         // 该光线与场景求交
//         Intersection obj2light = intersect(light);

//         // 如果反射击中光源
//         if (light2obj.happened && (light2obj.coords - lightPos).norm() < 1e-2)
//         {
//             // 获取改材质的brdf，这里的brdf为漫反射（brdf=Kd/pi）
//             Vector3f f_r = inter.m->eval(ray.direction, lightDir, N);
//             // 直接光照光 = 光源光 * brdf * 光线和物体角度衰减 * 光线和光源法线角度衰减 / 光线距离 / 该点的概率密度（1/该光源的面积）
//             L_dir = lightInter.emit * f_r * dotProduct(lightDir, N) * dotProduct(-lightDir, NN) / lightDistance / pdf_light;
//         }

//         // 俄罗斯轮盘赌，确定是否继续弹射光线
//         if (get_random_float() < RussianRoulette)
//         {
//             // 获取半平面上的随机弹射方向
//             Vector3f nextDir = inter.m->sample(ray.direction, N).normalized();
//             // 定义弹射光线
//             Ray nextRay(objPos, nextDir);
//             // 获取相交点
//             Intersection nextInter = intersect(nextRay);
//             // 如果有相交，且是与物体相交
//             if (nextInter.happened && !nextInter.m->hasEmission())
//             {
//                 // 该点间接光= 弹射点反射光 * brdf * 角度衰减 / pdf(认为该点四面八方都接收到了该方向的光强，为1/(2*pi)) / 俄罗斯轮盘赌值(强度矫正值)
//                 float pdf = inter.m->pdf(ray.direction, nextDir, N);
//                 Vector3f f_r = inter.m->eval(ray.direction, nextDir, N);
//                 L_indir = castRay(nextRay, depth + 1) * f_r * dotProduct(nextDir, N) / pdf / RussianRoulette;
//             }
//         }
//         // 最后返回直接光照和间接光照
//         return L_dir + L_indir;
//     }
//     // 如果光线与场景无交点
//     return Vector3f(0, 0, 0);
// }



Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    Intersection hit = this->intersect(ray);
    if (hit.happened) // 也可能是光源
    {
        if (hit.m->hasEmission()) // 光源
        {
            // 直接计算光源到像素的颜色
            // std::cout << itsec.emit << std::endl;
            if (depth == 0)
                return hit.m->getEmission(); // * itsec.m->eval(ray.direction, -1.0f * ray.direction, itsec.normal) / pow(itsec.distance, 2);
            else
                return Vector3f(0, 0, 0);
        }
        Vector3f p = hit.coords;
        Vector3f w_o = ray.direction;
        Vector3f N = hit.normal;

        // 直接光照
        Vector3f L_dir = Vector3f(0.0f, 0.0f, 0.0f);
        Intersection itsecLight;
        float pdf_light = 0.0f;
        sampleLight(itsecLight, pdf_light);

        // 灯光表面法线
        Vector3f NN = normalize(itsecLight.normal);
        // 物体点坐标
        auto &objPos = hit.coords;
        // 光源点坐标
        auto &lightPos = itsecLight.coords;

        // 光源采样
        Vector3f dir = lightPos - objPos;
        Vector3f dir2Light = normalize(dir);
        // float lightDistance = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;

        // 该光线与场景求交
        Intersection obj2light = intersect(Ray(p, dir2Light));

        if (obj2light.happened && obj2light.m->hasEmission()) // 应该击中位置相同的吧(obj2light.coords - lightPos).norm() < 1e-2
        {
            // printf("Hit Light");
            // Vector3f w_s = normalize(p - itsecLight.coords);
            // L_dir = itsecLight.emit * hit.m->eval(w_o, w_s, NN) * dotProduct(-1.0f * w_s, NN) * dotProduct(w_s, N) / pow(obj2light.distance, 2) / pdf_light;
            Vector3f f_r = hit.m->eval(w_o, dir2Light, N); // hit.m->eval(w_o, dir2Light, N); // hit.m->eval(w_o, w_s, N)
            // Vector3f light = obj2light.m->getEmission();   //  obj2light.m->getEmission()
            if (hit.m->m_type == REFLECT || hit.m->m_type == REFRACT) // 完全反射/折射到光源
            {
                if (1.0f - f_r.norm() < 0.002f)
                    L_dir = obj2light.m->getEmission();
            }
            else
            {  
                L_dir = obj2light.m->getEmission() * f_r * dotProduct(-dir2Light, NN) * dotProduct(dir2Light, N) / pow(obj2light.distance, 2) / pdf_light;
            }
        }

        // 间接光照
        if (get_random_float() > RussianRoulette)
        {
            return L_dir;
        }
        Vector3f L_indir = Vector3f(0.0f, 0.0f, 0.0f);

        Vector3f w_i = normalize(hit.m->sample(w_o, N)); // 在表面采样下一个光线（现在的出射光线）
        Ray nextRay = Ray(p, w_i);

        if (hit.m->m_type == REFLECT)
        {
            L_indir = castRay(nextRay, depth + 1) / RussianRoulette;
        }
        else if (hit.m->m_type == REFRACT)
        {
            L_indir = castRay(nextRay, depth + 1) / RussianRoulette;
        }
        else
        {
            Intersection inter = this->intersect(nextRay);
            if (inter.happened && !inter.m->hasEmission())
            {
                float pdf = hit.m->pdf(ray.direction, w_i, N); // 在表面采样的概率分布（均匀分布）
                // Vector3f q = inter.coords;
                Vector3f f_r = hit.m->eval(ray.direction, w_i, N); // 当前递归的点处的贡献程度（入对出）
                // std::cout << diffuse << std::endl;
                L_indir = castRay(nextRay, depth + 1) * f_r * dotProduct(w_i, N) / pdf / RussianRoulette;
            }
        }
        return L_dir + L_indir;
    }
    return this->backgroundColor;
} 