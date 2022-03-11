// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>

Vector3f color[1000000];
int rst::rasterizer::my_index(float a,float b)
{
    if(!hash.count({a,b}))hash[{a,b}] = my_idx++;
    return hash[{a,b}];
}

rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}


static bool insideTriangle(float x, float y, const Vector3f* _v)
{   
    Vector3f P(x,y,1.0f);
    Vector3f A = _v[0];
    Vector3f B = _v[1];
    Vector3f C = _v[2];
    Vector3f AB = B-A;
    Vector3f BC = C-B;
    Vector3f CA = A-C;
    Vector3f AP = P-A;
    Vector3f BP = P-B;
    Vector3f CP = P-C;

    float z1 = AB.cross(AP).z();
    float z2 = BC.cross(BP).z();
    float z3 = CA.cross(CP).z();
    
    return (z1>0&&z2>0&&z3>0)||(z1<0&&z2<0&&z3<0);
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();
    float xl = std::min({v[0].x(),v[1].x(),v[2].x()});
    float xh = std::max({v[0].x(),v[1].x(),v[2].x()});
    float yl = std::min({v[0].y(),v[1].y(),v[2].y()});
    float yh = std::max({v[0].y(),v[1].y(),v[2].y()});
    // xl*=2.0;xh*=2.0;yl*=2.0;yh*=2.0;
    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle
    for(int x=(int)xl;x<xh;++x)
    {
        for(int y=(int)yl;y<yh;++y)
        {
            bool flag=0;
            Vector3f temp=Vector3f::Zero();
            for(float i=0;i<1;i+=0.5)
                for(float j=0;j<1;j+=0.5)
                {
                    float nx = x+i;
                    float ny = y+j;
                    if(!insideTriangle(nx,ny,t.v))continue;
                    flag=1;
                    auto [alpha,beta,gamma] = computeBarycentric2D(nx,ny,t.v);
                    float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                    float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                    z_interpolated *= w_reciprocal;
                    int idx = my_index(nx,ny);
                    if(z_interpolated>=depth_buf[idx]) continue;
                    depth_buf[idx] = z_interpolated;
                    color[idx]= t.getColor()/4.0;
                }
            if(flag)
            {
                for(float i=0;i<1;i+=0.5)
                    for(float j=0;j<1;j+=0.5)
                    {
                        float nx = x+i;
                        float ny = y+j;
                        int idx = my_index(nx,ny);
                        temp += color[idx];
                    }
                set_pixel(Eigen::Vector3f(x,y,0),temp);
            }
            // if(!insideTriangle(x,y,t.v))continue;

            // auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
            // float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
            // float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
            // z_interpolated *= w_reciprocal;
            // long long buf_index = get_index(x,y);
            // if(z_interpolated>=depth_buf[buf_index])continue;
            // depth_buf[buf_index] = z_interpolated;
            // set_pixel(Eigen::Vector3f(x,y,1),t.getColor());
        }
    }

    // xl = std::min({v[0].x(),v[1].x(),v[2].x()});
    // xh = std::max({v[0].x(),v[1].x(),v[2].x()});
    // yl = std::min({v[0].y(),v[1].y(),v[2].y()});
    // yh = std::max({v[0].y(),v[1].y(),v[2].y()});
    // for(int x =(int)xl;x<xh;++x)
    // {
    //     for(int y=(int)yl;y<yh;++y)
    //     {
    //         if(!insideTriangle(x,y,t.v,1))continue;
    //         set_pixel(Vector3f(x,y,1),color[x*2][y*2]+color[x*2][y*2+1]+color[x*2+1][y*2]+color[x*2+1][y*2+1]);
    //     }
    // }

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    // w*=2;
    // h*=2;
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width+x;

    // return (height*2-1-y)*width*2 + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on