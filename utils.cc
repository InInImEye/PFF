/**
 * Taken from pfc folder for geometry and mesh creation
 * Sparsh Sinha
*/

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

// Need for geometry

#include <deal.II/base/geometry_info.h>

// Needed for C++ output

#include <iostream>
#include <fstream>

// Map
#include <map>

// Needed for math functions

#include <cmath>

// Rid of cumbersomeness

using namespace dealii;

/// Print mesh info
 
template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
 
  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;
 
    std::cout << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        std::cout << pair.first << '(' << pair.second << " times) ";
      }
    std::cout << std::endl;
  }
 
  std::ofstream out(filename);
  GridOut       grid_out;
  grid_out.write_vtu(triangulation, out);
  std::cout << " written to " << filename << std::endl << std::endl;
}

void create_2d_grid(
  const std::vector<Point<2>> &vertices,
  const std::vector<
	std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
	&               vertex_indices,
  Triangulation<2> &coarse_grid)
{
	std::vector<CellData<2>> cells(vertex_indices.size());
	for (unsigned int i = 0; i < cells.size(); ++i)
	{
	  for (unsigned int j = 0; j < vertex_indices[i].size(); ++j)
		cells[i].vertices[j] = vertex_indices[i][j];
	}

coarse_grid.create_triangulation(vertices, cells, SubCellData());
}

/// Create quarter notch // default manifold id 42 Origin aound 0,0 but not specified since manifold not assigned

void create_quarter_circular_notch(const double r,
	const double rsq,
	const types::manifold_id curve,
	Triangulation<2> &tria)
{
	const std::vector<Point<2>> vertices_c {{r, 0.},
									{r + rsq, 0.},
									{r / sqrt(2), r / sqrt(2)},
									{r + rsq, r + rsq},
									{0., r},
									{0., r + rsq}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_c = {{{1, 3, 0, 2}},
						 {{5, 4, 3, 2}}};
	
	create_2d_grid(vertices_c, cell_vertices_c, tria);
	
	for (auto &cell : tria.active_cell_iterators())
		for( auto &face : cell->face_iterators())
			{
				const Point<2> face_center = face->center();
				if (std::abs(face_center[0]) <= r && std::abs(face_center[1]) <= r)
					face->set_all_manifold_ids(curve);
			}
}

/// Create semi circular notch

void create_semi_circular_notch(const double r,
	const double rsq,
	const types::manifold_id curve,
	Triangulation<2> &tria)
{
	Triangulation<2> tria1, tria2;
	
	create_quarter_circular_notch( r, rsq, curve, tria1);
	create_quarter_circular_notch( r, rsq, curve, tria2);
	GridTools::rotate(numbers::PI / 2, tria2);
	GridGenerator::merge_triangulations(tria1, tria2, tria, 1.0e-12, true);//loses/ keeps manifold id // cpp does not allow middle values to be skipped like fn(a,,b) but allows right most to be skipped: fn(a,b) if def fn(a= ,b = ,c=) 
	// Ok manifold id gets copied but not set
}

/// Creating first mesh

void grid_1()
{
	
	// Circular coarse mesh
	Triangulation<2> tria_c;
	const double r = 1;
	const double rsq = 1;
	
	
	
	const Point<2> center(0., 0.);
	const PolarManifold<2> center_manifold(center);
	
	// Square coarse mesh
	Triangulation<2> tria_sq;
	const std::vector<Point<2>> vertices_sq {{0., - rsq},
									{r / sqrt(2), -rsq},
									{r + rsq, -rsq},
									{r + rsq + rsq, -rsq},
									{0., 0.},
									{r / sqrt(2), 0.},
									{r + rsq, 0.},
									{r + rsq + rsq, 0.},
									{r + rsq, r + rsq -  (r / sqrt(2))},
									{r + rsq + rsq, r + rsq -  (r / sqrt(2))},
									{r + rsq, r + rsq},
									{r + rsq + rsq, r + rsq}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq = {{{0, 1, 4, 5}},
						 {{1, 2, 5, 6}},
						 {{3, 7, 2, 6}},
						 {{7, 9, 6, 8}},
						 {{11, 10, 9, 8}}};
	create_2d_grid(vertices_sq, cell_vertices_sq, tria_sq);
	
	
	tria_c.set_manifold(42, center_manifold);
	
	Triangulation<2> tria_c2;
	//GridTools::shift(Point<2> (1., 0.5), tria_c2);
	// Merge the two meshes
	
	Triangulation<2> tria; 
	create_semi_circular_notch(r, rsq, 42, tria);
	/*
	for (auto &cell : tria.active_cell_iterators())
		for( auto &face : cell->face_iterators())
			{
				const Point<2> face_center = face->center();
				if (std::abs(face_center[0]- 0) <= r && std::abs(face_center[1]- 2.) <= r)
					face->set_all_manifold_ids(42);
				//std::cout << face->get_boundary_id() << std::endl;
			}
	//tria_c.set_manifold(42, center_manifold);
	*/
	tria.set_manifold(42, center_manifold); // here still need to set
	
	tria.refine_global(4);
	
	print_mesh_info(tria,"test.vtu");
}

/// Set manifolds 3d 2d /// Function overload

template <int dim=2>
 void multiple_manifold_set(const Triangulation<2> &triatmp,
	 Triangulation<2> &tria,
	 const std::vector<const Point<2>*> &Pts,
	 const std::vector<types::manifold_id> &man,
	 const unsigned int n_slices = 2,
	 const double height = 1)
{
    tria.copy_triangulation(triatmp);
    //unsigned int man[4] = {40, 41, 42, 43};
    //auto &it = Pts.at(2)[0];
    //std::cout << Pts.at(0)[0][1] << std::endl;
    //const PolarManifold<2> cm0(Pts.at(0)[0]);//, cm1(Pts[1]), cm2(Pts[2]), cm3(Pts[3]);
	if (Pts.size() == man.size())
	{
		for (unsigned int i = 0; i < Pts.size(); i++)
		{
			const PolarManifold<2> cm0(Pts.at(i)[0]);
			tria.set_manifold(man.at(i), cm0);
		}
	}
	else
		std::cout<< "ERROR:"
				 << std::endl 
				 << "Wrong Input size to multiple_manifold_set " 
				 << std::endl;
	/*  tria.set_manifold(40, cm0);
	  tria.set_manifold(41, cm1);
	  tria.set_manifold(42, cm2);
	  tria.set_manifold(43, cm3);*/
}

template <int dim=3>
 void multiple_manifold_set(const Triangulation<2> &triatmp,
	 Triangulation<3> &tria,
	 const std::vector<const Point<2>*> &Pts,
	 const std::vector<types::manifold_id> &man,
	 const unsigned int n_slices = 2,
	 const double height = 1)
{
	 
    const Tensor<1, 3> axis({0., 0., 1.});
    GridGenerator::extrude_triangulation(triatmp, n_slices, height, tria, true, {40, 41, 42, 43, numbers::flat_manifold_id});
    if (Pts.size() == man.size())
	{
		for (unsigned int i = 0; i < Pts.size(); i++)
		{
			const Point<3> Pz0(Pts.at(i)[0][0], Pts.at(i)[0][1], 0.);
			const CylindricalManifold<3> cym0(axis, Pz0);
			tria.set_manifold(man.at(i), cym0);
		}
	}
	else
		std::cout<< "!!ERROR:"
				 << std::endl 
				 << "Wrong Input size to multiple_manifold_set " 
				 << std::endl;
    
    /*
    
    
    tria.set_manifold(40, cym0);
    tria.set_manifold(41, cym1);
    tria.set_manifold(42, cym2);
    tria.set_manifold(43, cym3);

*/
}

 void multiple_curve_manifold_set_2d(
	 Triangulation<2> &tria,
	 const std::vector<const Point<2>*> &Pts,
	 const std::vector<types::manifold_id> &man)
{
   	if (Pts.size() == man.size())
	{
		for (unsigned int i = 0; i < Pts.size(); i++)
		{
			const PolarManifold<2> cm0(Pts.at(i)[0]);
			tria.set_manifold(man.at(i), cm0);
		}
	}
	else
		std::cout<< "ERROR:"
				 << std::endl 
				 << "Wrong Input size to multiple_manifold_set " 
				 << std::endl;
}


/// Create I bar geometry and mesh

template <int dim>
void grid_Ibar(Triangulation<dim> &tria, const unsigned int n_slices = 2, const double height = 1)
{
	Triangulation<2> triatmp, qcn0, qcn1, qcn2, qcn3, bar, bartmp;
	const double r = 10.;
	const double rsq = 20.;
	const Point<2> P0(0., r + rsq),
				   P1(2 * (r + rsq), r + rsq),
				   P2(0., (5 * rsq)-r),
				   P3(2 * (r + rsq), (5 * rsq)-r);
	
	
	/*
		2 -- 3
		|    |
		0 -- 1
	*/
	
	//0
	create_quarter_circular_notch(r, rsq, 40, qcn0);
	GridTools::rotate(-numbers::PI / 2, qcn0);
	GridTools::shift(P0, qcn0);
	//1
	create_quarter_circular_notch(r, rsq, 41, qcn1);
	GridTools::rotate(numbers::PI, qcn1);
	GridTools::shift(P1, qcn1);
	//2
	create_quarter_circular_notch(r, rsq, 42, qcn2);	
	GridTools::shift(P2, qcn2);
	//3
	create_quarter_circular_notch(r, rsq, 43, qcn3);
	GridTools::rotate(numbers::PI / 2, qcn3);
	GridTools::shift(P3, qcn3);
	
	//bar
	const std::vector<Point<2>> vertices_sq {{0., 0},
									{rsq, 0.},
									{0, rsq},
									{rsq, rsq}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq = {{{0, 1, 2, 3}}};
	create_2d_grid(vertices_sq, cell_vertices_sq, bartmp);
	GridGenerator::replicate_triangulation(bartmp, {2, 3}, bar);
	GridTools::shift(Point<2> (r, rsq + r), bar);
	//merge
	GridGenerator::merge_triangulations({&qcn0, &qcn1, &qcn2, &qcn3, &bar},
                                    triatmp,
                                    1.0e-12,
                                    true,
                                    false);

	//manifold
	
    
    multiple_manifold_set<dim>(triatmp, tria, {&P0, &P1, &P2, &P3}, {40, 41, 42, 43}, n_slices, height);
    /*
    const Tensor<1, 3> axis({0., 0., 1.});
    const Point<3> Pz0(P0(0), P0(1), 0.),
                   Pz1(P1(0), P1(1), 0.),
                   Pz2(P2(0), P2(1), 0.),
                   Pz3(P3(0), P3(1), 0.);
    const CylindricalManifold<3> cym0(axis, Pz0),
                                 cym1(axis, Pz1),
                                 cym2(axis, Pz2),
                                 cym3(axis, Pz3);
   
    
    
    tria.set_manifold(40, cym0);
    tria.set_manifold(41, cym1);
    tria.set_manifold(42, cym2);
    tria.set_manifold(43, cym3);

*/
  
}

/// Creates asymmetric double notches geometry
template <int dim>
void grid_asymm_notches(Triangulation<dim> &tria,
                        const unsigned int n_slices = 2,
                        const double height = 1)
{
	Triangulation<2> triatmp, triac1, triac2, bartmp, bartmp1, bar, bar3L, bar3R, bar5L, bar5R;
	const double r = 2.5;
	const double wid = 9.;
	const Point<2> P0(- wid , 2 * r), P1(wid , - 2 * r);
	
	create_semi_circular_notch(r, r, 40, triac1);
	GridTools::rotate(- numbers::PI / 2, triac1);
	GridTools::shift(Point<2> (-wid, 2* r), triac1);
	
	create_semi_circular_notch(r, r, 41, triac2);
	GridTools::rotate( numbers::PI / 2, triac2);
	GridTools::shift(Point<2> (wid, - 2* r), triac2);
	
	const std::vector<Point<2>> vertices_sq {{0., 0},
									{wid - (2 * r), 0.},
									{0, (2 * r)},
									{wid - (2 * r), (2 * r)}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq = {{{0, 1, 2, 3}}};
	create_2d_grid(vertices_sq, cell_vertices_sq, bartmp);
	GridGenerator::replicate_triangulation(bartmp, {2, 10}, bar);
	GridTools::shift(Point<2> ((2 * r) - wid, - 10 * r), bar);
	
	const std::vector<Point<2>> vertices_sq1 {{0., 0},
									{(2 * r), 0.},
									{0, (2 * r)},
									{(2 * r), (2 * r)}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq1 = {{{0, 1, 2, 3}}};
	create_2d_grid(vertices_sq1, cell_vertices_sq1, bartmp1);
	
	GridGenerator::replicate_triangulation(bartmp1, {1, 3}, bar3L);
	bar3R.copy_triangulation(bar3L);
	GridTools::shift(Point<2> (-wid, 4 * r), bar3L);
	GridTools::shift(Point<2> (wid - (2 * r), - 10 * r), bar3R);
	
	GridGenerator::replicate_triangulation(bartmp1, {1, 5}, bar5L);
	bar5R.copy_triangulation(bar5L);
	GridTools::shift(Point<2> (-wid, - 10 * r), bar5L);
	GridTools::shift(Point<2> (wid - (2 * r), 0.), bar5R);
	
	GridGenerator::merge_triangulations( {&triac1, &triac2, &bar, &bar3L, &bar5L, &bar3R, &bar5R},
		 triatmp,
		 1.0e-12,
		 true,
		 false);
	
	multiple_manifold_set<dim>(triatmp, tria, {&P0, &P1}, {40, 41}, n_slices, height);
	
	//print_mesh_info(tria,"results/asymm_notches.vtu");
}

/// Creates symmetric double notches geometry 1 and 2
template <int dim>
void grid_symm_notches_r(Triangulation<dim> &tria,
                         const double r = 5.,
                         const double rsq = 4.,
                         const unsigned int n_slices = 2,
                         const double height = 1)
{
	//const double rsq = .004;
	const double halflen = 25.;
	Triangulation<2> triac1, triac2, triatmp, bartmp, bar1, bar2;
	const Point<2> P0(- r - rsq , 0.), P1(r + rsq, 0.);
	
	create_semi_circular_notch(r, rsq, 40, triac1);
	GridTools::rotate(- numbers::PI / 2, triac1);
	GridTools::shift(Point<2> (-(rsq + r),0.), triac1);
	
	create_semi_circular_notch(r, rsq, 41, triac2);
	GridTools::rotate(numbers::PI / 2, triac2);
	GridTools::shift(Point<2> ((rsq + r),0.), triac2);
	
	const std::vector<Point<2>> vertices_sq {{0., 0},
									{r + rsq, 0.},
									{0, (halflen - r - rsq) / 2},
									{r + rsq, (halflen - r - rsq) / 2}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq = {{{0, 1, 2, 3}}};
	create_2d_grid(vertices_sq, cell_vertices_sq, bartmp);
	
	GridGenerator::replicate_triangulation(bartmp, {2, 2}, bar1);
	GridTools::shift(Point<2> (-(rsq + r),0.), bar1);
	
	bar2.copy_triangulation(bar1);
	
	GridTools::shift(Point<2> (0.,r + rsq), bar1);
	GridTools::shift(Point<2> (0.,- halflen), bar2);
	
	GridGenerator::merge_triangulations({&triac1, &triac2, &bar1, &bar2}, triatmp, 1.0e-12, true, false);
	
	multiple_manifold_set<dim>(triatmp, tria, {&P0, &P1}, {40, 41}, n_slices, height);
	
	//tria.refine_global(2);
	//print_mesh_info(tria,"symm_notches.vtu");
}

/// Creates a single notch grid geometry
/*
void grid_single_notch()
{
	Triangulation<2> tria;
	
	const double sd = .005;
	const double gp = .00001;
	const std::vector<Point<2>> vertices_sq {{0., 0.},
									{0., - rsq},
									{r / sqrt(2), -rsq},
									{r + rsq, -rsq},
									{r + rsq + rsq, -rsq},
									{0., 0.},
									{r / sqrt(2), 0.},
									{r + rsq, 0.},
									{r + rsq + rsq, 0.},
									{r + rsq, r + rsq -  (r / sqrt(2))},
									{r + rsq + rsq, r + rsq -  (r / sqrt(2))},
									{r + rsq, r + rsq},
									{r + rsq + rsq, r + rsq}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq = {{{0, 1, 4, 5}},
						 {{1, 2, 5, 6}},
						 {{3, 7, 2, 6}},
						 {{7, 9, 6, 8}},
						 {{11, 10, 9, 8}}};
	create_2d_grid(vertices_sq, cell_vertices_sq, tria);
	tria.refine_global(2);
	print_mesh_info(tria,"single_notch.vtu");
	
	
}*/

/// Creates 2 diagonally corner meshes
template <int dim>
void grid_2_corner_notches(Triangulation<dim> &tria,
                           const unsigned int n_slices = 2,
                           const double height = 1)
{
	Triangulation<2> triatmp, triac1, triac2, bar1, bar2;
	const double r = 2.5;
	const double in = 1.;
	const double thet = std::acos(in / r);
	const double len = 5.;
	const types::manifold_id cm0 = 40;
	const types::manifold_id cm1 = 41;
	const Point<2> Pt0(in, in), Pt1((2*len) - in, (2*len) - in);
	
	const std::vector<Point<2>> vertices_c {{in + r * std::sin(thet), 0.},
									{len, 0.},
									{r * cos(numbers::PI / 4) + in, r* sin(numbers::PI / 4) + in},
									{len, len},
									{0., in + r * std::sin(thet)},
									{0., len}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_c = {{{1, 3, 0, 2}},
						 {{5, 4, 3, 2}}};
	
	create_2d_grid(vertices_c, cell_vertices_c, triac1);
	triac2.copy_triangulation(triac1);
	triac1.set_all_manifold_ids(1);
	
	for (auto &cell : triac1.active_cell_iterators())
		for( auto &face : cell->face_iterators())
			{
				const Point<2> face_center = face->center();
				if (std::abs(face_center[0] - Pt0(0)) <= r && std::abs(face_center[1]- Pt0(1)) <= r)
					face->set_all_manifold_ids(cm0);
			}
	
	triac2.set_all_manifold_ids(1);
	for (auto &cell : triac2.active_cell_iterators())
		for( auto &face : cell->face_iterators())
			{
				const Point<2> face_center = face->center();
				if (std::abs(face_center[0] - Pt0(0)) <= r && std::abs(face_center[1]- Pt0(1)) <= r)
					face->set_all_manifold_ids(cm1);
			}
	GridTools::rotate(numbers::PI, triac2);
	GridTools::shift(Point<2> (2*len, 2*len), triac2);
	
	
	const std::vector<Point<2>> vertices_sq {{0., 0},
									{len, 0.},
									{0, len},
									{len, len}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq = {{{0, 1, 2, 3}}};
	create_2d_grid(vertices_sq, cell_vertices_sq, bar1);
	bar2.copy_triangulation(bar1);
	GridTools::shift(Point<2> (len, 0.), bar1);
	GridTools::shift(Point<2> (0., len), bar2);
	
	GridGenerator::merge_triangulations({&triac1, &bar1, &triac2, &bar2}, triatmp, 1.0e-12, true, false);
	
	multiple_manifold_set<dim>(triatmp, tria, {&Pt0, &Pt1}, {cm0, cm1}, n_slices, height);
	
	TransfiniteInterpolationManifold<dim> inner_manifold;
	inner_manifold.initialize(tria);
	tria.set_manifold (1, inner_manifold);
	
	//tria.refine_global(3);
	//print_mesh_info(tria,"two_corner_notches.vtu");
}

/// Creates a CT test geometry and mesh
template <int dim>
void grid_CT(Triangulation<dim> &tria,
             const unsigned int n_slices = 2,
             const double height = 1)
{
	const double rx = 25.4;
	const double ry = 13.97;
	
	const Point<2> P0(rx, ry);
	const Point<2> P1(rx, -ry);
	
	Triangulation<2> triatmp;
	
	
	GridIn<2> gridin;
	gridin.attach_triangulation(triatmp);
	std::ifstream f("inputFiles/CTtest1.msh");
	gridin.read_msh(f);
	
	triatmp.set_all_manifold_ids_on_boundary(11, 40);
	triatmp.set_all_manifold_ids_on_boundary(2, 40);
	
	triatmp.set_all_manifold_ids_on_boundary(3, 41);
	triatmp.set_all_manifold_ids_on_boundary(13, 41);
	
	multiple_manifold_set<dim>(triatmp, tria, {&P0, &P1}, {40, 41}, n_slices, height);
	
	//tria.refine_global(2);
	//print_mesh_info(tria,"CT.vtu");
}

/// Creates L bar

void grid_Lbar()
{
	Triangulation<2> tria;
	const double r = 0.1;
	const double rsq = 0.1;
	
	const std::vector<Point<2>> vertices_sq {{0., - rsq},
									{r / sqrt(2), -rsq},
									{r + rsq, -rsq},
									{r + rsq + rsq, -rsq},
									{0., 0.},
									{r / sqrt(2), 0.},
									{r + rsq, 0.},
									{r + rsq + rsq, 0.},
									{r + rsq, r + rsq -  (r / sqrt(2))},
									{r + rsq + rsq, r + rsq -  (r / sqrt(2))},
									{r + rsq, r + rsq},
									{r + rsq + rsq, r + rsq}};
	const std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
		cell_vertices_sq = {{{0, 1, 4, 5}},
						 {{1, 2, 5, 6}},
						 {{3, 7, 2, 6}},
						 {{7, 9, 6, 8}},
						 {{11, 10, 9, 8}}};
	create_2d_grid(vertices_sq, cell_vertices_sq, tria);
	//tria.refine_global(2);
	//print_mesh_info(tria,"Lbar.vtu");
}

/// Creates I bar in 3D
/*
void grid_3d_Ibar()
{
	Triangulation<2> tria;
	print_mesh_info(tria,"threeD_Ibar.vtu");
}
*/
/// Main

