#include "particles_system.hpp"


particles_system::particles_system() : Blocks(rect_vect(vect(-1.,-1.),vect(1.,1.)), vect(0.2, 0.2), 100)
{
  force_enabled = false;
  force_center = vect(0., 0.);

  // place particles in the domain
  double r=0.02;
  int N=100;
  for(int i=0; i<N; ++i)
  {
    int row=11;
    P.push_back(particle(vect((i%row*2.0+1.0)*r-1.0, (i/row*2.0+1.0)*r-1.0), vect(0.1, 0.0), 0.01, r, 10000.0, 1+i%2*0, rgb(1.0, i%2, 0.0)));
  }

  transfer_data(P, X,V);

  t=0.0;
  dt=0.0001;
  g=vect(0.0, -10.0);
}
particles_system::~particles_system()
{}
std::vector<particle> particles_system::GetParticles() const {
  return P;
}
void particles_system::AddEnvObj(env_object* env) {
  ENVOBJ.push_back(std::unique_ptr<env_object>(env));
}
void particles_system::status(std::ostream& out)
{
  out<<"Particles system"<<std::endl<<"Particles number = "<<P.size()<<std::endl;
}
void particles_system::step()
{
  for(std::size_t i=0; i<P.size(); ++i)
  {
    auto& p = P[i].p;
    auto& v = P[i].v;
    p += v * dt;
  }

  auto F = RHS();

  for(std::size_t i=0; i<P.size(); ++i)
  {
    auto& v = P[i].v;
    auto& f = F[i];
    v += (f + v * (-0.1)) * dt;
    const double limit = 100.;
    if (v.length() > limit) {
      v *= limit / v.length();
    }
  }

  Blocks.arrange(X);
  t += dt;
}
void particles_system::SetForce(vect center, bool enabled) {
    force_center = center;
    force_enabled = enabled;
}
void particles_system::SetForce(vect center) {
    force_center = center;
}
void particles_system::SetForce(bool enabled) {
    force_enabled = enabled;
}

void particles_system::transfer_data(const vector<particle>& P, vector<vect>& X, vector<vect>& V)
{
  size_t N=P.size();
  X.resize(N);
  V.resize(N);
  for(size_t i=0; i<N; ++i)
  {
    X[i]=P[i].p;
    V[i]=P[i].v;
  }
}

std::vector<vect> particles_system::RHS() const
{
  std::vector<vect> F(P.size());

  for(size_t i=0; i<F.size(); ++i)
  {
    auto& f = F[i];
    auto& p = P[i].p;
    auto& v = P[i].v;
    auto& part = P[i];
    // gravity
    f=g*part.m;
    // environment objects
    for(size_t k=0; k<ENVOBJ.size(); ++k)
    {
      auto& obj=ENVOBJ[k];
      f+=obj->F(p,v,part.r,part.sigma);
    }
    // point force
    if (force_enabled) {
        const double intensity = 0.1;
        const vect r = p - force_center; 
        f += r * (intensity / std::pow(r.length(), 3));
    }

    // dissipation
    f -= v * 0.1;
  }

  // pairwise interactions
  for(int bj=0; bj<Blocks.B.msize().j; ++bj)
  for(int bi=0; bi<Blocks.B.msize().i; ++bi)
  {
    mindex m(bi,bj);
    const std::vector<int>& b1=Blocks.B[m];
    for(size_t w1=0; w1<b1.size(); ++w1)
    {
      int p1 = b1[w1]; // first particle

      for(size_t k=0; k<Blocks.NEAR.size(); ++k)
      {
        mindex mnear = m + Blocks.NEAR[k];
        if(Blocks.B.valid(mnear))
        {
          const std::vector<int>& b2=Blocks.B[mnear];
          for(size_t w2=0; w2<b2.size(); ++w2)
          {
            int p2=b2[w2]; // second particle
            if(p1!=p2 && (P[p1].layers_mask & P[p2].layers_mask))
            {
              F[p1] += F12(P[p1].p, P[p1].v, P[p2].p, P[p2].v, 0.5*(P[p1].sigma+P[p2].sigma), P[p1].r+P[p2].r);
            }
          }
        }
      }
    }
  }

  for(size_t i=0; i<F.size(); ++i)
  {
    F[i]*=1.0/P[i].m;
  }
  return F;
}

vect F12(vect p1, vect /*v1*/, vect p2, vect /*v2*/, double sigma, double R)
{
  const double alpha=12.0;
  const double beta=6.0;
   
  double eps=sigma/(pow(2.0,alpha)-pow(2.0,beta));

  double r=p1.dist(p2);
  double F=r>R?0.0:eps*(pow(R/r, alpha)-pow(R/r, beta));
  return (p1-p2)*(F/r);
}
