//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
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
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    // shou xian wo men fen wei dir and indir 
    // ta men yi jing ba guang xian fa chu lai le 
    // mei yi ge cong xiang su dian gao chu lai de guang xian zhui zong ying gai fen wei dir and indir
    // su ji de jiu shi da dao obj ran hou fan she ran hou zai dao obj zai fan she etc..
    
    auto inter_P = intersect(ray);
    if(!inter_P.happened)return Vector3f();
    if(inter_P.m->hasEmission())return inter_P.m->getEmission();//zhe li shi de guang ke jian
    auto p =inter_P.coords;// gai yi xia
    auto N = inter_P.normal.normalized();
    auto w0 = ray.direction;
    Intersection inter;
    float pdf_light;
    sampleLight(inter,pdf_light);
    auto x = inter.coords;
    auto ws = x-p;
    ws = ws.normalized();
    auto NN = inter.normal.normalized();  
    auto emit = inter.emit;
    auto inter_temp = intersect(Ray(p,ws));
    float dist = (x-p).x*(x-p).x+(x-p).y*(x-p).y+(x-p).z*(x-p).z;
    Vector3f l_dir =Vector3f(0.0,0.0,0.0);
    // float N=0.0;//what's this
    if(inter_temp.distance-(double)(x-p).norm()>-EPSILON*10)
        l_dir = inter_P.m->eval(w0,ws,N)*emit*dotProduct(ws,N)*dotProduct(-ws,NN)/dist/(pdf_light+0.00001);

    Vector3f l_indir =Vector3f(0.0,0.0,0.0);
    float P_RR = get_random_float();
    if(P_RR<RussianRoulette)
    {
        auto wi = inter_P.m->sample(w0,N).normalized();//zen me zhi dao wo men xian zai de cai liao shi shen me 
        auto pdf_tt=inter_P.m->pdf(w0,wi,N);
        auto inter_tt=intersect(Ray(p,wi));
        if(inter_tt.happened&&!inter_tt.m->hasEmission())
        {
            l_indir = inter_P.m->eval(w0,wi,N)*castRay(Ray(inter_P.coords,wi),depth-1)*dotProduct(wi,N)/(pdf_tt+0.00001)/RussianRoulette;
        }
    }
    return l_dir+l_indir;
}