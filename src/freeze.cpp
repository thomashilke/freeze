#include <iostream>

#include <tfel/tfel.hpp>
#include <parameter/parameter.hpp>

class freeze {
public:
  freeze(const std::string& parameter_filename) {
    p.read_from_file(parameter_filename);

    particle_radius         = p.get_value<double>("particle-radius");
    alumina_diffusivity     = p.get_value<double>("alumina-diffusivity-coefficient");

    alumina_density = p.get_value<double>("alumina-density");
    electrolyte_density = p.get_value<double>("electrolyte-density");
    alumina_heat_capacity = p.get_value<double>("alumina-heat-capacity");
    electrolyte_heat_capacity = p.get_value<double>("liquid-electrolyte-heat-capacity");

    electrolyte_solid_liquid_low_t = p.get_value<double>("electrolyte-solid-liquid-low-t");
    electrolyte_solid_liquid_high_t = p.get_value<double>("electrolyte-solid-liquid-high-t");
    electrolyte_solid_liquid_latent_heat = p.get_value<double>("electrolyte-solid-liquid-latent-heat");
    
    electrolyte_diffusivity = p.get_value<double>("liquid-electrolyte-diffusivity-coefficient");
    t_inj                   = p.get_value<double>("alumina-initial-temperature");
    t_init                  = p.get_value<double>("electrolyte-initial-temperature");
    t_end                   = p.get_value<double>("t-end");
    r_max                   = p.get_value<double>("r-max");
    n                       = p.get_value<int>   ("space-subdivisions");
    m                       = p.get_value<int>   ("time-subdivisions");
  }

  void run() {
    heat_cartesian_coordinates();
    heat_spherical_coordinates();    
  }
  
  
private:
  parameter::collection p;
  
  double particle_radius;
  
  double alumina_diffusivity;
  double electrolyte_diffusivity;
  
  double alumina_density;
  double electrolyte_density;
  double alumina_heat_capacity;
  double electrolyte_heat_capacity;

  double electrolyte_solid_liquid_low_t;
  double electrolyte_solid_liquid_high_t;
  double electrolyte_solid_liquid_latent_heat;

  double t_inj;
  double t_init;
  double t_end;
  double r_max;

  std::size_t n;
  std::size_t m;

private:
  double initial_condition(const double* x) const {
    if (*x < particle_radius)
      return t_inj;
    else
      return t_init;
  }

  double diffusivity(const double* x) const {
    if (*x < particle_radius)
      return alumina_diffusivity;
    else
      return electrolyte_diffusivity;
  }

  double rho_cp(const double* x) const {
    if (*x < particle_radius)
      return alumina_density * alumina_heat_capacity;
    else
      return electrolyte_density * electrolyte_heat_capacity;
  }

  double beta(double u) const {
    const double
      u_low(electrolyte_solid_liquid_low_t / (electrolyte_density * electrolyte_heat_capacity)),
      u_high(u_low + electrolyte_solid_liquid_latent_heat);
    const double solid_liquid_delta_t(electrolyte_solid_liquid_high_t - electrolyte_solid_liquid_high_t);

    if (u < u_low)
      return u / (electrolyte_density * electrolyte_heat_capacity);
    else if (u < u_high)
      return electrolyte_solid_liquid_low_t
        + (electrolyte_solid_liquid_high_t - electrolyte_solid_liquid_low_t)
        / electrolyte_solid_liquid_latent_heat * (u - u_low);
    else
      return electrolyte_solid_liquid_high_t
        + (u - u_high)
        / (electrolyte_density * electrolyte_heat_capacity);
  }

  double beta_lipschitz_constant() const {
    return std::max(1.0/(alumina_density * alumina_heat_capacity),
                    1.0/(electrolyte_density * electrolyte_heat_capacity));
  }
  
  static double spherical_element(const double* x) {
    return (*x) * (*x);
  }

  void heat_cartesian_coordinates() {
    using namespace std::placeholders;

    const double delta_t(t_end / m);
  
    using cell_type = cell::edge;
    using mesh_type = fe_mesh<cell_type>;
    using fe_type = cell_type::fe::lagrange_p1;
    using fes_type = finite_element_space<fe_type>;
    using element_type = fes_type::element;

    mesh_type mesh(gen_segment_mesh(0.0, r_max, n));
    //submesh<cell_type> dmesh(mesh.get_boundary_submesh());
    submesh<cell_type> dmesh(mesh);

    fes_type fes(mesh, dmesh);

    element_type u_p(projector::lagrange<fe_type>(
      std::bind(&freeze::initial_condition, this, _1), fes));

    std::function<double(const double*)> D(std::bind(&freeze::diffusivity, this, _1));
    std::function<double(const double*)> R(std::bind(&freeze::rho_cp, this, _1));
    
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto u(a.get_trial_function());
      auto v(a.get_test_function());

      a += integrate<quad::edge::gauss3>((1 / delta_t) * make_expr(R) * u * v
                                         + make_expr(D) * d<1>(u) * d<1>(v), mesh);
    }

    linear_form<fes_type> f(fes);
    for (std::size_t k(0); k < m; ++k) {
      exporter::ascii<mesh_type>("solution_cc_" + std::to_string(k) + ".dat",
                                 to_mesh_vertex_data<fe_type>(u_p));

    
      auto v(f.get_test_function());
    
      f.clear();
      f += integrate<quad::edge::gauss3>((1 / delta_t) * make_expr(R)
                                         * make_expr<fe_type>(u_p) * v,
                                         mesh);

      u_p = a.solve(f);
    }
    exporter::ascii<mesh_type>("solution_cc_" + std::to_string(m) + ".dat",
                               to_mesh_vertex_data<fe_type>(u_p));
  }

  void heat_state_transition_cartesian_coordinates() {
    using namespace std::placeholders;

    using cell_type = cell::edge;
    using mesh_type = fe_mesh<cell_type>;
    using fe_type = cell_type::fe::lagrange_p1;
    using fes_type = finite_element_space<fe_type>;
    using element_type = fes_type::element;

    
    const double delta_t(t_end / m);
    const double mu(1.0/beta_lipschitz_constant());
    
    std::function<double(double)> B(std::bind(&freeze::beta, this, _1));
    std::function<double(const double*)> D(std::bind(&freeze::diffusivity, this, _1));
    std::function<double(const double*)> R(std::bind(&freeze::rho_cp, this, _1));    

    /*
     *  Mesh and finite element space
     */
    mesh_type mesh(gen_segment_mesh(0.0, r_max, n));
    submesh<cell_type> dmesh(mesh);
    fes_type fes(mesh, dmesh);


    // Initial condition for the enthalpy
    element_type u(projector::lagrange<fe_type>(
      std::bind(&freeze::initial_condition, this, _1), fes));

    // Initial condition for the temperature
    element_type t(u.transform(B));

    
    // Build the bilinear form
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto theta(a.get_trial_function());
      auto v(a.get_test_function());
      
      a += integrate<quad::edge::gauss3>(theta * v
                                         + delta_t / mu * d<1>(theta)*d<1>(v)
                                         , mesh);
    }

    // Iteration for each timestep
    linear_form<fes_type> f(fes);
    for (std::size_t k(0); k < m; ++k) {
      exporter::ascii<mesh_type>("solution_cc_" + std::to_string(k) + ".dat",
                                 to_mesh_vertex_data<fe_type>(u));

      {
        auto v(f.get_test_function());
    
        f.clear();
        f += integrate<quad::edge::gauss3>(make_expr<fe_type>(t) * v
                                           , mesh);
      }

      t = a.solve(f);

      auto tmp(t);
      tmp -= u.transform(B);
      tmp *= mu;
      u += tmp;
      //u = operator+<fe_type>(u, operator*<fe_type>(mu, (operator-<fe_type>(t, u.transform(B)))));
    }
    exporter::ascii<mesh_type>("solution_cc_" + std::to_string(m) + ".dat",
                               to_mesh_vertex_data<fe_type>(u));
  }

  void heat_spherical_coordinates() {
    using namespace std::placeholders;
  
    const double delta_t(t_end / m);
  
    using cell_type = cell::edge;
    using mesh_type = fe_mesh<cell_type>;
    using fe_type = cell_type::fe::lagrange_p1;
    using fes_type = finite_element_space<fe_type>;
    using element_type = fes_type::element;

    mesh_type mesh(gen_segment_mesh(0.0, r_max, n));
    //submesh<cell_type> dmesh(mesh.get_boundary_submesh());
    submesh<cell_type> dmesh(mesh);

    fes_type fes(mesh, dmesh);

    element_type u_p(
      projector::lagrange<fe_type>(
        std::bind(&freeze::initial_condition, this, _1),
        fes));
  
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto u(a.get_trial_function());
      auto v(a.get_test_function());

      std::function<double(const double*)> D(std::bind(&freeze::diffusivity, this, _1));
      a += integrate<quad::edge::gauss3>((1 / delta_t) * make_expr(spherical_element)
                                         * u * v
                                         + make_expr(spherical_element) * make_expr(D)
                                         * d<1>(u) * d<1>(v),
                                         mesh);
    }

    linear_form<fes_type> f(fes);
    for (std::size_t k(0); k < m; ++k) {
      exporter::ascii<mesh_type>("solution_sc_" + std::to_string(k) + ".dat",
                                 to_mesh_vertex_data<fe_type>(u_p));

    
      auto v(f.get_test_function());
    
      f.clear();
      f += integrate<quad::edge::gauss3>((1 / delta_t) * make_expr(spherical_element)
                                         * make_expr<fe_type>(u_p) * v,
                                         mesh);

      u_p = a.solve(f);
    }
    exporter::ascii<mesh_type>("solution_sc_" + std::to_string(m) + ".dat",
                               to_mesh_vertex_data<fe_type>(u_p));
  }
  
};


int main(int argc, char *argv[]) {
  try {
    if (argc != 2)
      std::cout << "Usage: " << argv[0] << "  <configuration file>" << std::endl;
    else {
      freeze f(argv[1]);
      f.run();
    }
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
