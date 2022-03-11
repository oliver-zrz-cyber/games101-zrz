#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   std::string splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}
// gao dong zen me gou jian de wai jia SAH
BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if(splitMethod=="NAIVE"){
        if (objects.size() == 1) {
            // Create leaf _BVHBuildNode_
            node->bounds = objects[0]->getBounds();
            node->object = objects[0];
            node->left = nullptr;
            node->right = nullptr;
            return node;
        }
        else if (objects.size() == 2) {
            node->left = recursiveBuild(std::vector{objects[0]});
            node->right = recursiveBuild(std::vector{objects[1]});

            node->bounds = Union(node->left->bounds, node->right->bounds);
            return node;
        }
        else {// da yu 2 wo men jiu xu yao hua fen
            Bounds3 centroidBounds;
            for (int i = 0; i < objects.size(); ++i)
                centroidBounds =
                    Union(centroidBounds, objects[i]->getBounds().Centroid());
            int dim = centroidBounds.maxExtent();
            switch (dim) {
            case 0:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().x <
                        f2->getBounds().Centroid().x;
                });
                break;
            case 1:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().y <
                        f2->getBounds().Centroid().y;
                });
                break;
            case 2:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().z <
                        f2->getBounds().Centroid().z;
                });
                break;
            }

            auto beginning = objects.begin();
            auto middling = objects.begin() + (objects.size() / 2);
            auto ending = objects.end();

            auto leftshapes = std::vector<Object*>(beginning, middling);
            auto rightshapes = std::vector<Object*>(middling, ending);

            assert(objects.size() == (leftshapes.size() + rightshapes.size()));

            node->left = recursiveBuild(leftshapes);
            node->right = recursiveBuild(rightshapes);
            node->bounds = Union(node->left->bounds, node->right->bounds);
        }

        return node;
    }
    else 
    {
        std::cout<<"SAH"<<std::endl;
        if (objects.size() == 1) {
            // Create leaf _BVHBuildNode_
            node->bounds = objects[0]->getBounds();
            node->object = objects[0];
            node->left = nullptr;
            node->right = nullptr;
            return node;
        }
        else if (objects.size() == 2) {
            node->left = recursiveBuild(std::vector{objects[0]});
            node->right = recursiveBuild(std::vector{objects[1]});
            node->bounds = Union(node->left->bounds, node->right->bounds);
            return node;
        }
        else{
            int bnumber=10;
            float mi = 1e9;
            int id = -1;
            Bounds3 centroidBounds;
            for (int i = 0; i < objects.size(); ++i)
                centroidBounds =
                    Union(centroidBounds, objects[i]->getBounds().Centroid());
            int dim = centroidBounds.maxExtent();
            switch (dim) {
            case 0:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().x <
                        f2->getBounds().Centroid().x;
                });
                break;
            case 1:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().y <
                        f2->getBounds().Centroid().y;
                });
                break;
            case 2:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().z <
                        f2->getBounds().Centroid().z;
                });
                break;
            }
            for(int i = 1;i<bnumber;++i)
            {
                auto beginning = objects.begin();
                auto middling = objects.begin() + (objects.size() *i/bnumber);
                auto ending = objects.end();

                auto leftshapes = std::vector<Object*>(beginning, middling);
                auto rightshapes = std::vector<Object*>(middling, ending);
                Bounds3 left,right;
                for(int k=0;k<leftshapes.size();++k)
                    left = Union(left,leftshapes[k]->getBounds().Centroid());
                for(int k=0;k<rightshapes.size();++k)right = Union(right,rightshapes[k]->getBounds().Centroid());
                float cost = 0.125+(left.SurfaceArea()*leftshapes.size()+right.SurfaceArea()*rightshapes.size())/centroidBounds.SurfaceArea();
                if(cost<mi)
                {
                    mi = cost,id = i;
                }
            }
            auto beginning = objects.begin();
            auto ending = objects.end();
            auto splitpoint = objects.begin() + (objects.size() *id/bnumber);
            auto leftshapes = std::vector<Object*>(beginning, splitpoint);
            auto rightshapes = std::vector<Object*>(splitpoint, ending);
            assert(objects.size() == (leftshapes.size() + rightshapes.size()));
            node->left = recursiveBuild(leftshapes);
            node->right = recursiveBuild(rightshapes);

            node->bounds = Union(node->left->bounds, node->right->bounds);

            return node;
        }
    }
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}
// shi bu shi xu yao yi zhi zhao dao ye zi jie dian 
Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const//an zhao na ge lei lai gao
{
    // TODO Traverse the BVH to find intersection
    auto t = node->bounds;
    Intersection ans;
    std::array<int,3>dirIsNeg;
    for(int i=0;i<3;++i)dirIsNeg[i] = ray.direction_inv[i]<0;
    if(!t.IntersectP(ray,ray.direction_inv,dirIsNeg))return ans;
    if(!node->left&&!node->right)// qiu jiao 
    {
        if(node->object!=nullptr)
            return node->object->getIntersection(ray);
    }
    else{
        auto t1 = getIntersection(node->left,ray);
        auto t2 = getIntersection(node->right,ray);
        if(t1.happened&&t2.happened)return t1.distance<t2.distance? t1:t2;
        else if(t1.happened)return t1;
        else return t2;
    }
}