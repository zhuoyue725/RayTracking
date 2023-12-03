#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    Scene scene(784, 784); // 784,784 512,512

    Material *pink = new Material(DIFFUSE, Vector3f(0.0f));
    pink->Kd = Vector3f(0.7f, 0.25f, 0.1f); // 0.20 , 0.020  0.015
    Material *purple = new Material(DIFFUSE, Vector3f(0.0f));
    purple->Kd = Vector3f(0.63f, 0.065f, 0.05f); // 0.20 , 0.020  0.015
    Material *white2 = new Material(DIFFUSE, Vector3f(0.0f));
    white2->Kd = Vector3f(1.0f, 1.0f, 1.0f); // 0.20 , 0.020  0.015
    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f); // 0.20 , 0.020  0.015
    Material *green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f); //0.04456 0.1432  0.0289
    Material *white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);// 0.23077 0.226 0.2164
    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(0.8f);  // 0.2069

    Material *totalRefract = new Material(REFRACT, Vector3f(0.0f));
    totalRefract->Ks = Vector3f(1.0f, 1.0f, 1.0f); // 0.23077 0.226 0.2164
    totalRefract->Kd = Vector3f(0.0f, 0.0f, 0.0f); // Vector3f(0.725f, 0.71f, 0.68f); // 0.23077 0.226 0.2164
    totalRefract->ior = 3.0f;

    Material *totalReflect = new Material(REFLECT, Vector3f(0.0f));
    totalReflect->Kd = Vector3f(0.0f, 0.0f, 0.0f); // 0.23077 0.226 0.2164
    totalReflect->Ks = Vector3f(1.0f, 1.0f, 1.0f); // Vector3f(0.725f, 0.71f, 0.68f); // 0.23077 0.226 0.2164
    totalReflect->ior = 1000.0f;

    Material *microfacet = new Material(MICROFACET, Vector3f(0.0f));
    microfacet->Kd = Vector3f(0.3f, 0.8f, 0.6f);   // 0.23077 0.226 0.2164
    microfacet->Ks = Vector3f(0.8f, 0.8f, 0.8f);   // Vector3f(0.725f, 0.71f, 0.68f); // 0.23077 0.226 0.2164
    microfacet->ior = 1.5f;
    microfacet->roughness = 0.3f;

    Sphere sphere1(Vector3f(400, 130, 200), 130, totalRefract); // left refract
    Sphere sphere2(Vector3f(160, 100, 350), 100, totalReflect); // right mirror
    Sphere sphere3(Vector3f(60, 40  , 80), 50, microfacet); // right refract

    MeshTriangle floor("../models/cornellbox/floor.obj", white);
    // MeshTriangle shortbox("../models/cornellbox/shortbox.obj", microfacet);
    // MeshTriangle tallbox("../models/cornellbox/tallbox.obj", microfacet);
    MeshTriangle left("../models/cornellbox/left.obj", pink);
    MeshTriangle right("../models/cornellbox/right.obj", purple);
    MeshTriangle light_("../models/cornellbox/light.obj", light);

    scene.Add(&floor);
    // scene.Add(&shortbox);
    // scene.Add(&tallbox);
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);
    scene.Add(&sphere1);
    scene.Add(&sphere2);
    scene.Add(&sphere3);

    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}