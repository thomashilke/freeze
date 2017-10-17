#include <iostream>

#include <tfel/tfel.hpp>
#include <parameter/parameter.hpp>

class freeze {
  enum class bc_type { dirichlet, neumann };
  enum class coordinates_type { spherical, cartesian };
  
public:
  freeze(const std::string& parameter_filename) {
    std::map<std::string, bc_type> bc_type_enum_map;
    bc_type_enum_map["dirichlet"] = bc_type::dirichlet;
    bc_type_enum_map["neumann"]   = bc_type::neumann;

    std::map<std::string, coordinates_type> coordinates_type_enum_map;
    coordinates_type_enum_map["spherical"] = coordinates_type::spherical;
    coordinates_type_enum_map["cartesian"] = coordinates_type::cartesian;
    
    
    p.read_from_file(parameter_filename); {
      particle_radius         = p.get_value<double>("particle-radius");

      electrolyte_density = p.get_value<double>("electrolyte-density");
      liquid_electrolyte_diffusivity = p.get_value<double>("liquid-electrolyte-diffusivity-coefficient");
      liquid_electrolyte_heat_capacity = p.get_value<double>("liquid-electrolyte-heat-capacity");
      solid_electrolyte_diffusivity = p.get_value<double>("liquid-electrolyte-diffusivity-coefficient");
      solid_electrolyte_heat_capacity = p.get_value<double>("liquid-electrolyte-heat-capacity");
    
      alumina_density = p.get_value<double>("alumina-density");
      alumina_heat_capacity = p.get_value<double>("alumina-heat-capacity");
      alumina_diffusivity     = p.get_value<double>("alumina-diffusivity-coefficient");

      electrolyte_sl_low_t = p.get_value<double>("electrolyte-sl-low-t");
      electrolyte_sl_high_t = p.get_value<double>("electrolyte-sl-high-t");
      electrolyte_sl_latent_heat = p.get_value<double>("electrolyte-sl-latent-heat");
    
      temp_inj                = p.get_value<double>("particle-initial-temperature");
      temp_init               = p.get_value<double>("electrolyte-initial-temperature");
      time_end                = p.get_value<double>("time-end");
      domain_size             = p.get_value<double>("domain-size");
      n                       = p.get_value<int>   ("space-subdivisions");
      m                       = p.get_value<int>   ("time-subdivisions");

      left_bc_type = p.get_enum_value<bc_type>("left-bc-type", bc_type_enum_map);
      right_bc_type = p.get_enum_value<bc_type>("right-bc-type", bc_type_enum_map);
      
      left_bc_value = p.get_value<double>("left-bc-value");
      right_bc_value = p.get_value<double>("right-bc-value");

      coordinates = p.get_enum_value<coordinates_type>("coordinates", coordinates_type_enum_map);

      output_prefix = p.get_value<std::string>("output-prefix");
      output_transient_solution = p.get_value<bool>("output-transient-solution");
      output_final_solution = p.get_value<bool>("output-final-solution");
      output_transition = p.get_value<bool>("output-transition");
      output_beta_function = p.get_value<bool>("output-beta-function");
      output_neumann_exact_solution = p.get_value<bool>("output-neumann-exact-solution");
    }
  }

  void run() {
    if (output_beta_function)
      show_beta_function();
    if (output_neumann_exact_solution)
      show_neumann_exact_solution();
    stefan();
  }

  void show_beta_function() {
    const double u_l(electrolyte_sl_low_t * solid_electrolyte_diffusivity
                     / (solid_electrolyte_heat_capacity * electrolyte_density));
    const double u_h(u_l + electrolyte_sl_latent_heat * electrolyte_density);

    
    // Length scale to get a enthalpy range to plot over
    double characteristic_enthalpy(1.0);
    if (electrolyte_sl_latent_heat * electrolyte_density > 0.0)
      characteristic_enthalpy = electrolyte_sl_latent_heat * electrolyte_density;

    // Min and max of the range
    const double u_min(u_l - characteristic_enthalpy);
    const double u_max(u_h + characteristic_enthalpy);

    // samples number
    const std::size_t n(500); 

    std::ofstream beta_file(output_prefix + "-beta.dat", std::ios::out);
    for (std::size_t i(0); i < n; ++i) {
      const double u(u_min + (static_cast<double>(i) / static_cast<double>(n - 1) * (u_max - u_min)));
      beta_file << u << " " << beta(u) << std::endl;
    }
  }

  void show_neumann_exact_solution() {
    const std::size_t n(500);

    std::ofstream neumann_es_file(output_prefix + "-neumann-exact-solution.dat", std::ios::out);
    const double t(time_end);

    for (std::size_t i(0); i < n; ++i) {
      const double x((static_cast<double>(i) / static_cast<double>(n - 1) * domain_size));
      neumann_es_file << x << " " << neumann_exact_solution(t, x) << std::endl;
    }
  }
  
private:
  parameter::collection p;

  // physical parameters
  double electrolyte_density;
  double liquid_electrolyte_diffusivity;
  double liquid_electrolyte_heat_capacity;
  double solid_electrolyte_diffusivity;
  double solid_electrolyte_heat_capacity;

  double electrolyte_sl_low_t;
  double electrolyte_sl_high_t;
  double electrolyte_sl_latent_heat;
  
  double alumina_density;
  double alumina_heat_capacity;
  double alumina_diffusivity;

  // model parameters
  double particle_radius;

  double temp_inj;
  double temp_init;
  
  double time_end;
  double domain_size;

  bc_type left_bc_type;
  bc_type right_bc_type;

  double left_bc_value;
  double right_bc_value;

  coordinates_type coordinates;
  
  // numerical parameters
  std::size_t n;
  std::size_t m;

  // output
  std::string output_prefix;
  bool output_transient_solution;
  bool output_final_solution;
  bool output_transition;
  bool output_beta_function;
  bool output_neumann_exact_solution;

private:
  using cell_type = cell::edge;
  using mesh_type = fe_mesh<cell_type>;
  using fe_type = cell_type::fe::lagrange_p1;
  using fes_type = finite_element_space<fe_type>;
  using element_type = fes_type::element;

  double neumann_exact_solution(const double t, const double x) {
    const double gamma(0.768955338463582);
    const double alpha(1.0);

    const double X(std::sqrt(2.0 * gamma * t));

    if (x > X) {
      return 0.0;
    } else {
      return - alpha * std::exp(gamma / 2.0) * std::sqrt(M_PI * gamma / 2.0)
        * (std::erf(std::sqrt(gamma / 2.0))
           - std::erf(x / std::sqrt(4.0 * t)));
    }
  }
  
  static double cartesian_volume_element(const double* x) {
    return 1.0;
  }

  static double spherical_volume_element(const double* x) {
    return *x * *x;
  }
  
  double initial_condition(const double* x) const {
    if (*x < particle_radius)
      return temp_inj;
    else
      return temp_init;
  }

  double beta(double u) const {
    const double u_l(electrolyte_sl_low_t * solid_electrolyte_diffusivity
                     / (solid_electrolyte_heat_capacity * electrolyte_density));
    const double u_h(u_l + electrolyte_sl_latent_heat * electrolyte_density);

    const double liquid_slope(liquid_electrolyte_diffusivity
                              / (liquid_electrolyte_heat_capacity
                                 * electrolyte_density));
    const double solid_slope(solid_electrolyte_diffusivity
                              / (solid_electrolyte_heat_capacity
                                 * electrolyte_density));

    const double mushy_slope((electrolyte_sl_high_t - electrolyte_sl_low_t)
                             / (electrolyte_sl_latent_heat * electrolyte_density));

    if (u < u_l)
      return electrolyte_sl_low_t + (u - u_l) * solid_slope;
    else if (u <= u_h)
      return electrolyte_sl_low_t + (u - u_l) * mushy_slope;
    else
      return electrolyte_sl_high_t + (u - u_h) * liquid_slope;
  }

  double inverse_beta(double temp) const {
    const double u_l(electrolyte_sl_low_t * solid_electrolyte_diffusivity
                     / (solid_electrolyte_heat_capacity * electrolyte_density));
    const double u_h(u_l + electrolyte_sl_latent_heat * electrolyte_density);

    const double liquid_slope(liquid_electrolyte_diffusivity
                              / (liquid_electrolyte_heat_capacity
                                 * electrolyte_density));
    const double solid_slope(solid_electrolyte_diffusivity
                              / (solid_electrolyte_heat_capacity
                                 * electrolyte_density));

    double mushy_slope((electrolyte_sl_high_t - electrolyte_sl_low_t)
                       / (electrolyte_sl_latent_heat * electrolyte_density));
    if (mushy_slope < 1e-6)
      mushy_slope = 1.0;
    
    if (temp < electrolyte_sl_low_t)
      return u_l + (temp - electrolyte_sl_low_t) / solid_slope;
    else if (temp < electrolyte_sl_high_t)
      return u_l + (temp - electrolyte_sl_low_t) / mushy_slope;
    else
      return u_h + (temp - electrolyte_sl_high_t) / liquid_slope;
  }

  double beta_lipschitz_constant() const {
    const double liquid_slope(liquid_electrolyte_diffusivity
                              / (liquid_electrolyte_heat_capacity
                                 * electrolyte_density));
    const double solid_slope(solid_electrolyte_diffusivity
                              / (solid_electrolyte_heat_capacity
                                 * electrolyte_density));

    const double mushy_slope((electrolyte_sl_high_t - electrolyte_sl_low_t)
                             / (electrolyte_sl_latent_heat * electrolyte_density));
    
    return std::max({solid_slope, mushy_slope, liquid_slope});
  }
  
  static double spherical_element(const double* x) {
    return (*x) * (*x);
  }

  void setup_boundary_conditions(fes_type& fes) {
    const mesh_type& mesh(fes.get_mesh());
    submesh<cell_type> dm(mesh.get_boundary_submesh());

    if (left_bc_type == bc_type::dirichlet
        and right_bc_type == bc_type::dirichlet) {

      fes.set_dirichlet_boundary_condition(dm, [this](const double* x){
          return *x < this->domain_size ? this->left_bc_value : this->right_bc_value;
        });
    } else if (left_bc_type == bc_type::dirichlet) {
      submesh<cell_type> left_b(dm.query_cells([this](const double* x){ return *x < this->domain_size / 2.0; }));
      fes.set_dirichlet_boundary_condition(left_b, [this](const double* x){ return this->left_bc_value; });
    } else if (right_bc_type == bc_type::dirichlet) {
      submesh<cell_type> right_b(dm.query_cells([this](const double* x){ return *x > this->domain_size / 2.0; }));
      fes.set_dirichlet_boundary_condition(right_b, [this](const double* x){ return this->right_bc_value; });
    }
  }

  std::set<double> build_level_set(const mesh_data<double, fe_mesh<cell_type> >& data, double level) {
    std::set<double> level_crossings;

    const auto& mesh(data.get_mesh());
    const auto& vertices(mesh.get_vertices());
    const auto& cells(mesh.get_cells());
    const auto& values(data.get_values());
    
    for (std::size_t i(0); i < mesh.get_cell_number(); ++i) {
      const double v1(values.at(cells.at(i, 0), 0)), v2(values.at(cells.at(i, 1), 0));
      const double x1(vertices.at(cells.at(i, 0), 0)), x2(vertices.at(cells.at(i, 1), 0));
      
      if ((v1 > level and v2 < level)
          or (v1 < level and v2 > level)) {
        level_crossings.insert(x1 + (x2 - x1) / (v2 - v1) * (level - v1));
      } else if (v1 == level
                 and (v2 > level or v2 < level)) {
        level_crossings.insert(x1);
      } else if (v2 == level
                 and (v1 > level or v1 < level)) {
        level_crossings.insert(x2);
      } else if (v1 == level and v2 == level) {
        level_crossings.insert(x1);
        level_crossings.insert(x2);
      } else if ((v1 > level and v2 > level)
                 or (v1 < level and v2 < level)) {
        // no transition in this element, do nothing
      } else {
        throw parameter::string_builder("build_level_set: This shohldn't happen, i missed a case. (level = ")
          (level)(", v1 = ")(v1)(", v2 = ")(v2)(")").str();
      }
    }

    return level_crossings;
  }

  
  void stefan() {
    using namespace std::placeholders;

    const double delta_t(time_end / m);
    const double mu(1.0/beta_lipschitz_constant());

    std::function<double(double)> B(std::bind(&freeze::beta, this, _1));
    std::function<double(double)> Binv(std::bind(&freeze::inverse_beta, this, _1));

    std::function<double(const double*)> volume_element;
    switch (coordinates) {
    case coordinates_type::cartesian:
      volume_element = cartesian_volume_element;
      break;
    case coordinates_type::spherical:
      volume_element = spherical_volume_element;
      break;
    }


    /*
     *  Mesh and finite element space
     */
    mesh_type mesh(gen_segment_mesh(0.0, domain_size, n));
    submesh<cell_type> empty_boundary(mesh);
    submesh<cell_type> dm(mesh.get_boundary_submesh());
    fes_type fes(mesh, empty_boundary);

    setup_boundary_conditions(fes);


    // Initial condition for the enthalpy
    element_type temp(projector::lagrange<fe_type>(
      std::bind(&freeze::initial_condition, this, _1), fes));

    
    // Initial condition for the temperature
    element_type u(temp.transformed(Binv));

    
    // Build the bilinear form
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto theta(a.get_trial_function());
      auto v(a.get_test_function());
      
      a += integrate<quad::edge::gauss3>(make_expr(volume_element) * theta * v
                                         + (delta_t / mu) * make_expr(volume_element) * d<1>(theta) * d<1>(v)
                                         , mesh);

    }

    // Storage for the transition location at each timestep
    std::vector<std::pair<double,double> > transitions;
    
    // Iteration for each timestep
    linear_form<fes_type> f(fes);
    for (std::size_t k(0); k < m; ++k) {
      if (output_transient_solution) {
        exporter::ascii<mesh_type>(output_prefix + "-ts-" + std::to_string(k) + ".dat",
                                   to_mesh_vertex_data<fe_type>(temp));
      }

      /*
       *  Build the rhs
       */
      {
        auto v(f.get_test_function());
        
        f.clear();
        f += integrate<quad::edge::gauss3>(make_expr(volume_element) * make_expr<fe_type>(temp) * v
                                           , mesh);
        /*
         *  Handle the Neumann boundary conditions
         */
        if (left_bc_type == bc_type::neumann) {
          submesh<cell_type> left_b(dm.query_cells([this](const double* x){ return *x < this->domain_size / 2.0; }));
          f += integrate<quad::point::eval>(- left_bc_value * delta_t / mu * make_expr(volume_element) * v, left_b);
        }

        if (right_bc_type == bc_type::neumann) {
          submesh<cell_type> right_b(dm.query_cells([this](const double* x){ return *x > this->domain_size / 2.0; }));
          f += integrate<quad::point::eval>(right_bc_value * delta_t / mu * make_expr(volume_element) * v, right_b);
        }
      }


      /*
       *  Chernoff formula
       */
      const auto theta(a.solve(f));

      auto tmp(theta);
      tmp -= temp;
      tmp *= mu;
      u += tmp;
      temp = u.transformed(B);

      /*
       * Post processing to get the interface
       */
      const double transition_temp((electrolyte_sl_low_t + electrolyte_sl_high_t) / 2.0);
      std::set<double> transition(build_level_set(to_mesh_vertex_data<fe_type>(temp), transition_temp));

      if (transition.size())
        transitions.push_back(std::make_pair(delta_t * k, *transition.begin()));
    }
    
    if (output_final_solution) {
      exporter::ascii<mesh_type>(output_prefix + ".dat",
                                 to_mesh_vertex_data<fe_type>(temp));

    }

    if (output_transition) {
      std::ofstream transition_file(output_prefix + "-transitions.dat",
                                    std::ios::out);
      for (const auto& t: transitions)
        transition_file << t.first << " " << t.second << std::endl;
    }
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
