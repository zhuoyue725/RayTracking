//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include <mutex>
#include <thread>
#include <pthread.h>
#include "Scene.hpp"
#include "Renderer.hpp"


inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }
// void saveFrameBuffer(std::vector<Vector3f> &framebuffer, int w, int h);
void saveFrameBuffer(Vector3f *framebuffer, int w, int h);
// void renderPixels(const Scene &scene, int start, int end, float scale, float imageAspectRatio, int spp, std::vector<Vector3f> &framebuffer, Vector3f eye_pos);
void renderPixels(const Scene &scene, int start, int end, float scale, float imageAspectRatio, int spp, Vector3f* framebuffer, Vector3f eye_pos , const int t);
// void renderPixels(const Scene &scene, int start, int end, float scale, float imageAspectRatio, int spp, std::vector<Vector3f> &framebuffer, Vector3f &eye_pos);
void saveToBuffer(std::vector<Vector3f> &framebuffer, std::vector<Vector3f> &colorsBuffer, int start, int startIdx, int end);
void saveToFile(Vector3f *framebuffer, const Scene &scene);

const float EPSILON = 0.00001;
int finishPixels = 0;
float invWhole = 0.0f;
// create a mutex to protect framebuffer
// std::mutex mtx;

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    // std::vector<Vector3f> framebuffer(scene.width * scene.height); // 动态内存，多线程访问容易冲突
    Vector3f *framebuffer = new Vector3f[scene.width * scene.height];

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;
    invWhole = 1.0f / (float)(scene.height *scene.width);
    // change the spp value to change sample ammount for each pixel
    int spp = 256;
    std::cout << "SPP: " << spp << "\n";

    // create a vector to store thread objects
    std::vector<std::thread> threads;

    // decide how many threads to use and how many pixels per thread
    int numThreads = 12; // std::thread::hardware_concurrency();          // get the number of available cores
    int pixelsPerThread = scene.width * scene.height / numThreads; // divide the total pixels by the number of threads
    std::cout << "numThreads: " << numThreads << "\n";
    std::cout << "maxNumThreads: " << std::thread::hardware_concurrency() << "\n";
    std::cout << "pixelsPerThread: " << pixelsPerThread << "\n";

    // create threads and assign subtasks to them
    for (int t = 0; t < numThreads; ++t)
    {
        // calculate the start and end index for each subtask
        int start = t * pixelsPerThread;
        int end = (t == numThreads - 1) ? scene.width * scene.height : (t + 1) * pixelsPerThread;
        // create a new thread with renderPixels function and pass the scene and range as arguments
        threads.push_back(std::thread(renderPixels, std::ref(scene), start, end, scale, imageAspectRatio, spp, framebuffer, eye_pos, t));
    }
    // 主线程实时保存数据
    saveToFile(framebuffer,scene);
    // wait for all threads to finish
    for (auto &thread : threads)
    // for (int t = 0; t < numThreads; ++t)
    {
            thread.join();
    }

    // for (uint32_t j = 0; j < scene.height; ++j) 
    // {
    //     for (uint32_t i = 0; i < scene.width; ++i) 
    //     {
    //         // generate primary ray direction
    //         float x = (2 * (i + 0.5) / (float)scene.width - 1) *
    //                   imageAspectRatio * scale;
    //         float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

    //         Vector3f dir = normalize(Vector3f(-x, y, 1));
    //         for (int k = 0; k < spp; k++)
    //         {
    //             framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;  
    //         }
    //         m++;
    //     }
    //     UpdateProgress(j / (float)scene.height);
    //     // std::cout << (j % 5) << std::endl;
    //     if ((int)(j % 5) == 0)
    //     {        
    //         // std::thread t(saveFrameBuffer, framebuffer,scene.width,scene.height);
    //         saveFrameBuffer(framebuffer,scene.width,scene.height);
    //         // 线程后台完成
    //         // t.detach();
    //     }
    // }

    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
void saveFrameBuffer(Vector3f *framebuffer, int w, int h)
{
    // save framebuffer to file
    FILE *fp = fopen("binary.ppm", "wb");
    // printf("Saving");
    (void)fprintf(fp, "P6\n%d %d\n255\n", w, h);
    for (auto i = 0; i < h * w; ++i)
    {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);
    // printf("Save Finish");
}
// define a function that takes a range of pixels and computes their colors
void renderPixels(const Scene &scene, int start, int end, float scale, float imageAspectRatio, int spp, Vector3f* framebuffer, Vector3f eye_pos,const int t)
{
    // printf("Thread %d to %d\n",start,end);
    // int total = end - start;

    int workStartHeight = start / scene.width;
    // std::vector<Vector3f> colorsBuffer(total);
    // int lastidx = 0, lastm = start;
    for (int m = start; m < end; ++m)
    {
        // get the pixel coordinates from the index
        int i = m % scene.width;
        int j = m / scene.width;
        // generate primary ray direction
        float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                  imageAspectRatio * scale;
        float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

        Vector3f dir = normalize(Vector3f(-x, y, 1));
        for (int k = 0; k < spp; k++)
        {
            framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
        }
        finishPixels += 1;
        // std::cout << (j % 5) << std::endl;
        
        // if ((int)((j - workStartHeight+1) % 10) == 0) // 做了5行
        // {
        //     // printf("Saving m = %d j = %d mod = %d\n", m,j,(j - workStartHeight + 1));
        //     // std::thread t(saveFrameBuffer, framebuffer,scene.width,scene.height);
        //     saveFrameBuffer(framebuffer, scene.width, scene.height); // 保存文件
        // }
    }
}
void saveToFile(Vector3f *framebuffer, const Scene &scene)
{
    float process = finishPixels * invWhole; // 进度
    while (process < 0.99f)
    {
        process = finishPixels * invWhole; // 进度
        UpdateProgress(process);
        if ((int)process % 5 == 0)
        {
            saveFrameBuffer(framebuffer, scene.width, scene.height); // 只有主线程可以保存文件
        }
        std::this_thread::sleep_for(std::chrono::seconds(3)); // 等待3秒
    }
}