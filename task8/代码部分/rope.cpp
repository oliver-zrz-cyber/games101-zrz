#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        bool flag=false;
        Mass* pre_ptr=nullptr;
        for(int i=0;i<num_nodes;++i)
        {
            double t = (double)i/(double)(num_nodes-1);
            Mass *ptr = new Mass(start+(end-start)*t,node_mass,false);
            masses.push_back(ptr);
        }
        for(int i=1;i<num_nodes;++i)
        {
            Spring *ptr = new Spring(masses[i-1],masses[i],k);
            springs.push_back(ptr);
        }
        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
       }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            auto a = s->m1;
            auto b = s->m2;
            s->m1->forces +=s->k*(((b->position-a->position).norm()-s->rest_length)*(b->position-a->position)/(b->position-a->position).norm());
            s->m2->forces -=s->k*(((b->position-a->position).norm()-s->rest_length)*(b->position-a->position)/(b->position-a->position).norm());
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces+=gravity*m->mass;
                m->velocity = m->velocity+m->forces/m->mass*delta_t;
                m->position = m->position + m->velocity*delta_t;
                // TODO (Part 2): Add global damping
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            if(s->k==-1)continue;
            auto a = s->m1;
            auto b = s->m2;
            s->m1->forces +=s->k*(((b->position-a->position).norm()-s->rest_length)*(b->position-a->position)/(b->position-a->position).norm());
            s->m2->forces -=s->k*(((b->position-a->position).norm()-s->rest_length)*(b->position-a->position)/(b->position-a->position).norm());
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                m->forces+=gravity*m->mass;
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                auto a = m->forces/m->mass;
                float dumping = 0.00005;

                m->position = m->position+ (1-dumping)*(m->position - m->last_position)+a*delta_t*delta_t;
                m->last_position = temp_position;
                // TODO (Part 4): Add global Verlet damping
            }
            m->forces = Vector2D(0,0);
        }
        for(auto &s:springs)
        {
            if(s->k!=-1)continue;
            auto d = s->m2->position-s->m1->position;
            auto d_norm = d.norm();
            auto d_dir = d/d_norm;
            auto offset1 = d_dir*0.5*(d_norm-s->rest_length);
            auto offset2 = -d_dir*0.5*(d_norm-s->rest_length);

            if(s->m1->pinned)offset2*=2,offset1 = Vector2D(0,0);
            if(s->m2->pinned)offset1*=2,offset2 = Vector2D(0,0);

            s->m1->position += offset1;
            s->m2->position += offset2;

        }
    }
}
