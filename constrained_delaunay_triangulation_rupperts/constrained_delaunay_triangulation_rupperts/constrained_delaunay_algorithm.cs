using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace constrained_delaunay_triangulation
{
    public static class constrained_delaunay_algorithm
    {
        // A Delaunay Refinement Algorithm for Quality 2-Dimensional Mesh Generation
        // Jim Ruppert (NASA Ames Research Center)
        // https://www.cis.upenn.edu/~cis610/ruppert.pdf
        // Delaunay Refinement Algorithms for Triangular Mesh Generation
        // Jonathan Richard Shewchuk jrs@cs.berkeley.edu (May 21, 2001)
        // https://people.eecs.berkeley.edu/~jrs/papers/2dj.pdf
        // Guranteed-Quality Mesh Generation for curved surfaces L.Paul Chew Cornell University Ithaca, NY
        // https://kogs-www.informatik.uni-hamburg.de/~tchernia/SR_papers/chew93.pdf
        // https://ecommons.cornell.edu/handle/1813/6899

        // There Are Planar Graphs Almost as Good as the Complete Graph L.PAUL CHEW
        // https://core.ac.uk/download/pdf/82457584.pdf


        // Q- Morph Algorithm References
        // ] Steven J. Owen, Matthew L. Staten, Scott A. Canann, and Sunil Saigal.
        // Advancing Front Quadrilateral Meshing Using Triangle Transformations.In Proceedings of the 7th International Meshing Roundtable.October 1998.
        // http://www.imr.sandia.gov/papers/imr7/owen98.ps.gz
        // https://www.researchgate.net/profile/Steven_Owen3/publication/221561698_Advancing_Front_Quadrilateral_Meshing_Using_Triangle_Transformations/links/0046352cee8df0718b000000/Advancing-Front-Quadrilateral-Meshing-Using-Triangle-Transformations.pdf
        // https://pdfs.semanticscholar.org/6984/9f35237f416811e23a776167da7a6aabaeab.pdf
        // http://pages.cs.wisc.edu/~csverma/CS899_09/qmorph.pdf



        public static double eps = 0.000001; // 10^-6
        private static mesh_store main_mesh;

        public static void create_constrained_mesh(int the_surface_index, List<int> inner_surface_indices, ref List<pslg_datastructure.surface_store> the_surface_data)
        {
            // step : 1 Create outter and inner edges of the surface
            List<pslg_datastructure.edge2d> outter_edges = new List<pslg_datastructure.edge2d>(); // variable to store outter edges
            bool seed_outter_edge_exists = false;
            // set the outter encapsulating nodes;
            outter_edges = set_surface_edges(the_surface_data[the_surface_index], ref seed_outter_edge_exists);

            // set the inner encapsulating nodes
            List<pslg_datastructure.edge2d>[] inner_edges = new List<pslg_datastructure.edge2d>[inner_surface_indices.Count]; // variable to store inner surfaces edges
            bool[] seed_inner_edge_exists = new bool[inner_surface_indices.Count];
            for (int i = 0; i < inner_surface_indices.Count; i++)
            {
                // set the inner encapsulating nodes one at a time;
                inner_edges[i] = new List<pslg_datastructure.edge2d>();
                seed_inner_edge_exists[i] = false;
                inner_edges[i] = set_surface_edges(the_surface_data[inner_surface_indices[i]], ref seed_inner_edge_exists[i]);
            }
            // ________________________________________________________________________________________________________________________________

            // step : 2 Add the edges vertex to point list
            List<pslg_datastructure.point2d> all_pts = new List<pslg_datastructure.point2d>(); // variable to store all the node
            foreach (pslg_datastructure.edge2d o_ed in outter_edges) // add the outter edges pts
            {
                all_pts.Add(new pslg_datastructure.point2d(all_pts.Count, o_ed.end_pt.x, o_ed.end_pt.y)); // just add the end point of the edges (Since the edges are closed => adding end point alone covers all the surface nodes)
            }
            for (int i = 0; i < inner_surface_indices.Count; i++) // add the inner edges pts
            {
                foreach (pslg_datastructure.edge2d i_ed in inner_edges[i]) // add the i_th inner edges pts
                {
                    all_pts.Add(new pslg_datastructure.point2d(all_pts.Count, i_ed.end_pt.x, i_ed.end_pt.y)); // just add the end point of the edges (Since the edges are closed => adding end point alone covers all the surface nodes)
                }
            }
            // ________________________________________________________________________________________________________________________________

            bool is_encroched = false;
            do
            {
                // step : 3 Encroched subsegment (Encroached subsegments are given priority over skinny triangles)
                // Refinement operation 1 => split a segment by adding a vertex at its midpoint
                is_encroched = false;
                // step : 3A Check whether the outter edge is encroched
                List<pslg_datastructure.point2d> amend_nodes = new List<pslg_datastructure.point2d>();
                encroched_segment_recursive_divide(ref outter_edges, outter_edges.Count - 1, 0, ref all_pts, ref amend_nodes, ref is_encroched, ref seed_outter_edge_exists);
                // step : 3B Check whether the inner edges are encroched
                for (int i = 0; i < inner_surface_indices.Count; i++) // cycle through all the inner surfaces
                {
                    encroched_segment_recursive_divide(ref inner_edges[i], inner_edges[i].Count - 1, 0, ref all_pts, ref amend_nodes, ref is_encroched, ref seed_inner_edge_exists[i]);
                }
                // ________________________________________________________________________________________________________________________________
                all_pts.AddRange(amend_nodes);

            } while (is_encroched == true);

            // step : 4 Compute initial Delaunay triangulation
            // delaunay triangulation using Bowyer watson incremental algorithm
            main_mesh = new mesh_store();
            main_mesh.Add_multiple_points(all_pts);
            main_mesh.Finalize_mesh(the_surface_data[the_surface_index]);
            // ________________________________________________________________________________________________________________________________

            // step : xx Well sized triangle condition is mean d
            List<pslg_datastructure.triangle2d> temp_sized_triangle = main_mesh.local_output_triangle.OrderBy(obj => obj.circum_radius).ToList();
            double mean_size = temp_sized_triangle[0].circum_radius * 20;


            do
            {
                // step: 5 Grade the triangle with circum radius to shortest edge ratio and find the worst triangle
                mesh_store.triangle_store worst_triangle = get_worst_triangle(main_mesh.all_triangles, the_surface_data[the_surface_index], outter_edges, inner_edges);

                // step : 6 Refine the mesh with circum radius to shortest edge ratio
                // and To produce an-edge of length less than h, the radius of circumcircle must be less than h.
                if (worst_triangle != null && (worst_triangle.circumradius_shortest_edge_ratio > Form1.the_static_class.B_var || worst_triangle.circum_radius > mean_size))
                {
                    pslg_datastructure.point2d inner_surface_pt = new pslg_datastructure.point2d(-1, worst_triangle.circum_center.x, worst_triangle.circum_center.y);
                    // step : 6A Refine the outter edges which are encroched by the failed triangles circum center
                    is_encroched = false;
                    encroched_segment_single_divide(ref outter_edges, ref inner_surface_pt, ref is_encroched, ref seed_outter_edge_exists);

                    if (is_encroched == true)
                    {
                        // incremental add point
                        main_mesh.Add_single_point(inner_surface_pt);
                        goto loopend;
                    }

                    // step : 6B Refine all the inner edges which are encroched by the failed triangles circum center
                    for (int j = 0; j < inner_surface_indices.Count; j++)
                    {
                        // cycle through all the inner surfaces
                        is_encroched = false;
                        encroched_segment_single_divide(ref inner_edges[j], ref inner_surface_pt, ref is_encroched, ref seed_inner_edge_exists[j]);

                        if (is_encroched == true)
                        {
                            // incremental add point
                            main_mesh.Add_single_point(inner_surface_pt);
                            goto loopend;
                        }
                    }

                    // Refinement operation 2 => split a triangle by adding a vertex at its circum center
                    // step : 6C Refine all the inner surface with the failed triangles circum center
                    main_mesh.Add_single_point(inner_surface_pt);

                loopend:;
                    // Finalize the mesh (to remove the faces out of bounds)
                    // main_mesh.Finalize_mesh(the_surface_data[the_surface_index]);

                }
                else
                {
                    // Exit if no bad triangles
                    break;
                }

            } while (true);

            // Finalize the mesh (to remove the faces out of bounds)
            main_mesh.Finalize_mesh(the_surface_data[the_surface_index]);

            // ________________________________________________________________________________________________________________________________

            // step : 7 Add the mesh to the surface
            the_surface_data[the_surface_index].my_mesh = new pslg_datastructure.mesh2d();
            the_surface_data[the_surface_index].is_meshed = true;
            the_surface_data[the_surface_index].my_mesh = new pslg_datastructure.mesh2d(main_mesh.local_input_points, main_mesh.local_output_edges, main_mesh.local_output_triangle);
            // ________________________________________________________________________________________________________________________________

            // step : 8 set the mesh seeds to the surface edges
            // step : 8A set the mesh seeds to the outter surface edges
            if (seed_outter_edge_exists == false)
            {
                the_surface_data[the_surface_index].encapsulating_seed_edges = new List<pslg_datastructure.edge2d>();
                the_surface_data[the_surface_index].encapsulating_seed_edges.AddRange(outter_edges);
            }

            // step : 8B set the mesh seeds to the inner surface edges
            for (int i = 0; i < inner_surface_indices.Count; i++) // cycle through all the inner surfaces
            {
                if (seed_inner_edge_exists[i] == false)
                {
                    the_surface_data[inner_surface_indices[i]].encapsulating_seed_edges = new List<pslg_datastructure.edge2d>();
                    the_surface_data[inner_surface_indices[i]].encapsulating_seed_edges.AddRange(inner_edges[i]);
                }
            }
            // ________________________________________________________________________________________________________________________________
            // End
        }

        public static List<pslg_datastructure.edge2d> set_surface_edges(pslg_datastructure.surface_store the_surface, ref bool seed_control)
        {
            List<pslg_datastructure.edge2d> the_edges = new List<pslg_datastructure.edge2d>(); // the edges is the return variable
            if (the_surface.encapsulating_seed_edges.Count == 0)
            {
                // if the encapsulating seed nodes are not found then the outter surface is not meshed. So set the seed nodes as surface boundary nodes
                seed_control = false; // seed doesn't exists
                foreach (pslg_datastructure.edge2d ed in the_surface.surface_edges)
                {
                    the_edges.Add(ed); // surface boundary nodes are added 
                }
            }
            else
            {
                seed_control = true; // seed exists
                foreach (pslg_datastructure.edge2d ed in the_surface.encapsulating_seed_edges)
                {
                    the_edges.Add(ed); // seed nodes are found which means the outter surface is already meshed and use the nodes of that mesh
                }
            }
            return the_edges;
        }

        public static void encroched_segment_recursive_divide(ref List<pslg_datastructure.edge2d> surface_edges, int end_index, int start_index, ref List<pslg_datastructure.point2d> all_nodes, ref List<pslg_datastructure.point2d> amend_nodes, ref bool to_continue, ref bool seed_exists)
        {
            if (seed_exists == false)
            {
                // Recorsion algorithm to split the segements
                for (int i = end_index; i >= start_index; i--) // go reverse to avoid overflow error
                {
                    // Form a diametral circle
                    pslg_datastructure.point2d cicle_center = new pslg_datastructure.point2d(-1, surface_edges[i].mid_pt.x, surface_edges[i].mid_pt.y); // circle center
                    double circle_radius_squared = Math.Pow((surface_edges[i].edge_length * 0.5) - eps, 2);


                    foreach (pslg_datastructure.point2d pt in all_nodes)
                    {


                        // (x - center_x)^2 + (y - center_y)^2 < radius^2 (then the point is inside the circle)
                        if ((Math.Pow((pt.x - cicle_center.x), 2) + Math.Pow((pt.y - cicle_center.y), 2)) < circle_radius_squared)
                        {
                            // create two segements from one segment by splitting at the mid point
                            pslg_datastructure.edge2d seg_1 = new pslg_datastructure.edge2d(surface_edges[i].edge_id, surface_edges[i].start_pt, surface_edges[i].mid_pt); // create a segment 1 with start_pt to mid_pt
                            pslg_datastructure.edge2d seg_2 = new pslg_datastructure.edge2d(surface_edges.Count, surface_edges[i].mid_pt, surface_edges[i].end_pt); // create a segment 2 with mid_pt to end_pt
                            amend_nodes.Add(new pslg_datastructure.point2d(all_nodes.Count + amend_nodes.Count, surface_edges[i].mid_pt.x, surface_edges[i].mid_pt.y));

                            surface_edges.RemoveAt(i); // remove the edge which encroches

                            // Insert the newly created two element edge list to the main list (at the removed index)
                            surface_edges.Insert(i, seg_1);
                            surface_edges.Insert(i + 1, seg_2);
                            to_continue = true;

                            // Recursion occurs here
                            encroched_segment_recursive_divide(ref surface_edges, i + 1, i, ref all_nodes, ref amend_nodes, ref to_continue, ref seed_exists);

                            break; // break is very important to avoid no longer continuing the search for the pt encroching i_th index segment
                        }
                    }
                }
            }

        }

        public static void encroched_segment_single_divide(ref List<pslg_datastructure.edge2d> surface_edges, ref pslg_datastructure.point2d circumcenter_pt, ref bool to_continue, ref bool seed_exists)
        {
            if (seed_exists == false)
            {
                // algorithm to split the segement
                for (int i = surface_edges.Count - 1; i >= 0; i--) // go reverse to avoid overflow error
                {
                    // Form a diametral circle
                    pslg_datastructure.point2d cicle_center = new pslg_datastructure.point2d(-1, surface_edges[i].mid_pt.x, surface_edges[i].mid_pt.y); // circle center
                    double circle_radius_squared = Math.Pow((surface_edges[i].edge_length * 0.5) - eps, 2);


                    // (x - center_x)^2 + (y - center_y)^2 < radius^2 (then the point is inside the circle)
                    if ((Math.Pow((circumcenter_pt.x - cicle_center.x), 2) + Math.Pow((circumcenter_pt.y - cicle_center.y), 2)) < circle_radius_squared)
                    {
                        // create two segements from one segment by splitting at the mid point
                        pslg_datastructure.edge2d seg_1 = new pslg_datastructure.edge2d(surface_edges[i].edge_id, surface_edges[i].start_pt, surface_edges[i].mid_pt); // create a segment 1 with start_pt to mid_pt
                        pslg_datastructure.edge2d seg_2 = new pslg_datastructure.edge2d(surface_edges.Count, surface_edges[i].mid_pt, surface_edges[i].end_pt); // create a segment 2 with mid_pt to end_pt
                        circumcenter_pt = new pslg_datastructure.point2d(-1, surface_edges[i].mid_pt.x, surface_edges[i].mid_pt.y);

                        surface_edges.RemoveAt(i); // remove the edge which encroches

                        // Insert the newly created two element edge list to the main list (at the removed index)
                        surface_edges.Insert(i, seg_1);
                        surface_edges.Insert(i + 1, seg_2);
                        to_continue = true;

                        break; // break is very important to avoid no longer continuing the search for the pt encroching i_th index segment
                    }
                }
            }
        }

        public static mesh_store.triangle_store get_worst_triangle(List<mesh_store.triangle_store> all_triangles, pslg_datastructure.surface_store the_surface, List<pslg_datastructure.edge2d> outter_edges, List<pslg_datastructure.edge2d>[] inner_edges)
        {
            // get the worst triangle
            // condition is the circume center must lies inside the surface or inside the diametrical circle

            List<mesh_store.triangle_store> graded_triangle = all_triangles.OrderBy(obj => obj.circumradius_shortest_edge_ratio).ThenBy(obj => obj.circum_radius).ToList();
            mesh_store.triangle_store worst_triangle = null;

            // step: 5B select the worst triangle by thether the triangle's circum center lies inside the surface 
            for (int i = graded_triangle.Count - 1; i >= 0; i--)
            {
                if (is_point_encroched(outter_edges, new pslg_datastructure.point2d(-1, graded_triangle[i].circum_center.x, graded_triangle[i].circum_center.y)) == true)
                {
                    worst_triangle = graded_triangle[i];
                    goto loopend;
                }
                else if (the_surface.pointinsurface(graded_triangle[i].circum_center.x, graded_triangle[i].circum_center.y) == true)
                {
                    worst_triangle = graded_triangle[i];
                    goto loopend;
                }
                else
                {
                    for (int j = 0; j < inner_edges.Length; j++)
                    {
                        if (is_point_encroched(inner_edges[j], new pslg_datastructure.point2d(-1, graded_triangle[i].circum_center.x, graded_triangle[i].circum_center.y)) == true)
                        {
                            worst_triangle = graded_triangle[i];
                            goto loopend;
                        }
                    }
                }
            }
        loopend:;
            return worst_triangle;
        }

        public static bool is_point_encroched(List<pslg_datastructure.edge2d> surface_edges, pslg_datastructure.point2d circumcenter_pt)
        {
            for (int i = surface_edges.Count - 1; i >= 0; i--) // go reverse to avoid overflow error
            {
                // Form a diametral circle
                pslg_datastructure.point2d cicle_center = new pslg_datastructure.point2d(-1, surface_edges[i].mid_pt.x, surface_edges[i].mid_pt.y); // circle center
                double circle_radius_squared = Math.Pow((surface_edges[i].edge_length * 0.5) - eps, 2);


                // (x - center_x)^2 + (y - center_y)^2 < radius^2 (then the point is inside the circle)
                if ((Math.Pow((circumcenter_pt.x - cicle_center.x), 2) + Math.Pow((circumcenter_pt.y - cicle_center.y), 2)) < circle_radius_squared)
                {
                    return true;
                }
            }

            return false;
        }

        public static void delete_constrained_mesh(int the_surface_index, List<int> inner_surface_indices, ref List<pslg_datastructure.surface_store> the_surface_data)
        {
            // delete the mesh
            the_surface_data[the_surface_index].is_meshed = false;
            the_surface_data[the_surface_index].my_mesh = new pslg_datastructure.mesh2d();

            // find the outter mesh surface index using the surface id
            int this_surface_encapsulating_surface_id = the_surface_data[the_surface_index].encapsulating_surface_id;
            int encapsulating_surface_index = the_surface_data.FindIndex(obj => obj.surface_id == this_surface_encapsulating_surface_id);

            //  check if the outter mesh surface index exists
            if (encapsulating_surface_index != -1)
            {
                if (the_surface_data[encapsulating_surface_index].is_meshed == false)
                {
                    the_surface_data[the_surface_index].encapsulating_seed_edges = new List<pslg_datastructure.edge2d>();
                }
            }

            // delete the inner mesh seed (only seeds)
            for (int i = 0; i < inner_surface_indices.Count; i++) // cycle through all the inner surfaces
            {
                if (the_surface_data[inner_surface_indices[i]].is_meshed == false)
                {
                    the_surface_data[inner_surface_indices[i]].encapsulating_seed_edges = new List<pslg_datastructure.edge2d>();
                }
            }

        }

        #region " Delaunay Triangulation - Bowyer Watson Incremental Algorithm"
        public class mesh_store
        {
            //#### Output Variables #######
            public List<pslg_datastructure.point2d> local_input_points = new List<pslg_datastructure.point2d>();
            public List<pslg_datastructure.edge2d> local_output_edges = new List<pslg_datastructure.edge2d>();
            public List<pslg_datastructure.triangle2d> local_output_triangle = new List<pslg_datastructure.triangle2d>();
            //#####################################

            // Local Variables
            private List<point_store> _all_points = new List<point_store>();
            private List<edge_store> _all_edges = new List<edge_store>();
            private List<triangle_store> _all_triangles = new List<triangle_store>();
            public bool is_meshed = false;

            private List<int> unique_edgeid_list = new List<int>();
            private List<int> unique_triangleid_list = new List<int>();

            // super triangle points
            private point_store s_p1, s_p2, s_p3;

            public List<point_store> all_points
            {
                get { return this._all_points; }
            }

            public List<triangle_store> all_triangles
            {
                get { return this._all_triangles; }
            }

            public List<edge_store> all_edges
            {
                get { return this._all_edges; }
            }

            public mesh_store()
            {
                // Empty constructor
            }

            public void Add_multiple_points(List<pslg_datastructure.point2d> parent_inpt_points)
            {
                // Call this first
                this._all_points = new List<point_store>();
                // Transfer the parent variable to local point_store list
                foreach (pslg_datastructure.point2d pt in parent_inpt_points)
                {
                    this._all_points.Add(new point_store(pt.id, pt.x, pt.y, pt));
                }

                // set the points to local list
                local_input_points = parent_inpt_points;

                // sort the points first by x and then by y
                this._all_points = this._all_points.OrderBy(obj => obj.x).ThenBy(obj => obj.y).ToList();


                // intialize the edges and triangles
                this._all_edges = new List<edge_store>();
                this._all_triangles = new List<triangle_store>();

                // Create an imaginary triangle that encloses all the point set
                set_bounding_triangle(this.all_points);

                foreach (point_store i_pt in this.all_points)
                {
                    // incemental add point
                    incremental_point_addition(i_pt);
                }
            }

            public void Add_single_point(pslg_datastructure.point2d parent_onemore_point)
            {
                // dont call this before calling Add_multiple_points
                // Add the point to local list
                this._all_points.Add(new point_store(this.all_points.Count, parent_onemore_point.x, parent_onemore_point.y, parent_onemore_point));
                // Add to local point list
                this.local_input_points.Add(this.all_points[this.all_points.Count - 1].get_parent_data_type);

                // call the incremental add point
                incremental_point_addition(this.all_points[this.all_points.Count - 1]);
            }

            private void incremental_point_addition(point_store pt)
            {
                // Find the index of triangle containing this point
                int containing_triangle_index = this.all_triangles.FindIndex(obj => obj.is_point_inside(pt) == true);

                if (containing_triangle_index != -1)
                {
                    // Point is inside a triangle
                    point_store tri_a = this.all_triangles[containing_triangle_index].pt1;
                    point_store tri_b = this.all_triangles[containing_triangle_index].pt2;
                    point_store tri_c = this.all_triangles[containing_triangle_index].pt3;

                    // remove the single triangle
                    remove_triangle(containing_triangle_index);

                    // add the three triangles
                    int first_triangle_id, second_triangle_id, third_triangle_id;

                    first_triangle_id = add_triangle(pt, tri_a, tri_b);
                    second_triangle_id = add_triangle(pt, tri_b, tri_c);
                    third_triangle_id = add_triangle(pt, tri_c, tri_a);

                    // Flip the bad triangles recursively
                    flip_bad_edges(first_triangle_id, pt);
                    flip_bad_edges(second_triangle_id, pt);
                    flip_bad_edges(third_triangle_id, pt);
                }
                else
                {
                    // Point lies on the edge
                    // Find the edge which is closest to the pt
                    int the_edge_id = this.all_edges.OrderBy(obj => obj.find_closest_distance(pt)).ToArray()[0].edge_id;
                    int the_edge_index = this.all_edges.FindIndex(obj => obj.edge_id == the_edge_id);


                    int first_tri_index = this.all_triangles.FindLastIndex(obj => obj.contains_edge(this.all_edges[the_edge_index]));
                    int second_tri_index = this.all_triangles.FindIndex(obj => obj.contains_edge(this.all_edges[the_edge_index]));

                    // the other point of first and second triangle sharing the common edge
                    point_store first_tri_other_pt = this.all_triangles[first_tri_index].get_other_pt(this.all_edges[the_edge_index]);
                    point_store second_tri_other_pt = this.all_triangles[second_tri_index].get_other_pt(this.all_edges[the_edge_index]);
                    // start and end point of the edge
                    point_store edge_a_pt = this.all_edges[the_edge_index].start_pt;
                    point_store edge_b_pt = this.all_edges[the_edge_index].end_pt;

                    // Remove the common edge
                    unique_edgeid_list.Add(this.all_edges[the_edge_index].edge_id);
                    this._all_edges.RemoveAt(the_edge_index);

                    // remove the two triangle
                    remove_triangle(first_tri_index);
                    remove_triangle(second_tri_index);

                    // add the four triangles
                    int first_triangle_id, second_triangle_id, third_triangle_id, fourth_triangle_id;

                    first_triangle_id = add_triangle(pt, edge_a_pt, first_tri_other_pt);
                    second_triangle_id = add_triangle(pt, first_tri_other_pt, edge_b_pt);
                    third_triangle_id = add_triangle(pt, edge_b_pt, second_tri_other_pt);
                    fourth_triangle_id = add_triangle(pt, second_tri_other_pt, edge_a_pt);

                    // Flip the bad triangles recursively
                    flip_bad_edges(first_triangle_id, pt);
                    flip_bad_edges(second_triangle_id, pt);
                    flip_bad_edges(third_triangle_id, pt);
                    flip_bad_edges(fourth_triangle_id, pt);
                }
            }

            public void Finalize_mesh(pslg_datastructure.surface_store the_surface)
            {
                // Call this after calling Add_multiple_points & Add_single_point (if required)
                // Finalize mesh is the final step to transfer the mesh from local data to global
                local_output_edges = new List<pslg_datastructure.edge2d>();
                local_output_triangle = new List<pslg_datastructure.triangle2d>();


                int i = 0;
                foreach (triangle_store the_tri in this.all_triangles)
                {
                    // Check if the triangles lies inside the surface
                    pslg_datastructure.triangle2d temp_tri = new pslg_datastructure.triangle2d(i, the_tri.pt1.get_parent_data_type, the_tri.pt2.get_parent_data_type, the_tri.pt3.get_parent_data_type);

                    if (the_surface.pointinsurface(temp_tri.shrunk_vertices[0].x, temp_tri.shrunk_vertices[0].y) == false ||
                        the_surface.pointinsurface(temp_tri.shrunk_vertices[1].x, temp_tri.shrunk_vertices[1].y) == false ||
                        the_surface.pointinsurface(temp_tri.shrunk_vertices[2].x, temp_tri.shrunk_vertices[2].y) == false)
                    {
                        // continue because the face is not in surface
                        continue;
                    }
                    local_output_triangle.Add(temp_tri);
                    i++;
                }

                i = 0;
                foreach (edge_store the_edge in this.all_edges)
                {
                    // Check if the edges lies inside the surface
                    pslg_datastructure.edge2d temp_edge = new pslg_datastructure.edge2d(i, the_edge.start_pt.get_parent_data_type, the_edge.end_pt.get_parent_data_type);

                    if (local_output_triangle.Exists(obj => obj.edge_exists(temp_edge) == true)) //so check whether this edge is associated with a face
                    {
                        local_output_edges.Add(temp_edge);
                        i++;
                    }
                }

                // set the mesh complete
                is_meshed = true;
            }


            private void remove_triangle(int tri_index)
            {
                int edge_index1 = this.all_edges.FindLastIndex(obj => obj.edge_id == this.all_triangles[tri_index].e1.edge_id);
                int edge_index2 = this.all_edges.FindLastIndex(obj => obj.edge_id == this.all_triangles[tri_index].e2.edge_id);
                int edge_index3 = this.all_edges.FindLastIndex(obj => obj.edge_id == this.all_triangles[tri_index].e3.edge_id);

                // Edge 1
                if (edge_index1 != -1)
                {
                    this._all_edges[edge_index1].remove_triangle(this.all_triangles[tri_index]);
                }
                // Edge 2
                if (edge_index2 != -1)
                {
                    this._all_edges[edge_index2].remove_triangle(this.all_triangles[tri_index]);
                }
                // Edge 3
                if (edge_index3 != -1)
                {
                    this._all_edges[edge_index3].remove_triangle(this.all_triangles[tri_index]);
                }

                unique_triangleid_list.Add(this.all_triangles[tri_index].tri_id);
                this._all_triangles.RemoveAt(tri_index);
            }

            private int add_triangle(point_store p1, point_store p2, point_store p3)
            {
                int edge_index1 = -1;
                int edge_index2 = -1;
                int edge_index3 = -1;

                //if (p1.Equals(p2) || p2.Equals(p3) || p3.Equals(p1))
                //{
                //    int k = 10;
                //    k = 100;
                //}

                // Edge 1
                edge_index1 = this.all_edges.FindLastIndex(obj => obj.Equals_without_orientation(new edge_store(-1, p1, p2)));
                if (edge_index1 == -1)
                {
                    edge_index1 = this.all_edges.Count;
                    this._all_edges.Add(new edge_store(get_unique_edge_id(), p1, p2));
                }

                // Edge 2
                edge_index2 = this.all_edges.FindLastIndex(obj => obj.Equals_without_orientation(new edge_store(-1, p2, p3)));
                if (edge_index2 == -1)
                {
                    edge_index2 = this.all_edges.Count;
                    this._all_edges.Add(new edge_store(get_unique_edge_id(), p2, p3));
                }

                // Edge 3
                edge_index3 = this.all_edges.FindLastIndex(obj => obj.Equals_without_orientation(new edge_store(-1, p3, p1)));
                if (edge_index3 == -1)
                {
                    edge_index3 = this.all_edges.Count;
                    this._all_edges.Add(new edge_store(get_unique_edge_id(), p3, p1));
                }

                // Triangle
                this._all_triangles.Add(new triangle_store(get_unique_triangle_id(), p1, p2, p3, this.all_edges[edge_index1], this.all_edges[edge_index2], this.all_edges[edge_index3]));

                // Update the triangle details to the edge
                this._all_edges[edge_index1].add_triangle(this.all_triangles[this.all_triangles.Count - 1]);
                this._all_edges[edge_index2].add_triangle(this.all_triangles[this.all_triangles.Count - 1]);
                this._all_edges[edge_index3].add_triangle(this.all_triangles[this.all_triangles.Count - 1]);

                return this.all_triangles[this.all_triangles.Count - 1].tri_id;
            }

            private void flip_bad_edges(int tri_id, point_store pt)
            {
                // find the triangle with input id
                int tri_index = this.all_triangles.FindLastIndex(obj => obj.tri_id == tri_id);
                // find the edge of this triangle whihc does not contain pt
                int common_edge_index = this.all_edges.FindLastIndex(obj => obj.edge_id == this.all_triangles[tri_index].get_other_edge_id(pt));

                if (common_edge_index == -1)
                {
                    return;
                }

                // main method to legalize edges by recursively flipping all the illegal edges
                int neighbour_tri_index = this.all_triangles.FindIndex(obj => obj.tri_id == this.all_edges[common_edge_index].other_triangle_id(this.all_triangles[tri_index]));
                //  int neighbour_tri_index = common_edge_index != -1 ? this.all_triangles.FindIndex(obj => obj.tri_id == this.all_edges[common_edge_index].other_triangle_id(this.all_triangles[tri_index])) : -1;

                // legalize only if the triangle has a neighbour
                if (neighbour_tri_index != -1)
                {
                    // check whether the newly added pt is inside the neighbour triangle
                    if (this.all_triangles[neighbour_tri_index].is_point_in_circumcircle(pt) == true)
                    {
                        // find the other vertex of the closest triangle
                        point_store neighbour_tri_other_vertex = this.all_triangles[neighbour_tri_index].get_other_pt(this.all_edges[common_edge_index]);
                        point_store edge_a_pt = this.all_edges[common_edge_index].start_pt;
                        point_store edge_b_pt = this.all_edges[common_edge_index].end_pt;

                        // Remove the common edge
                        unique_edgeid_list.Add(this.all_edges[common_edge_index].edge_id);
                        this._all_edges.RemoveAt(common_edge_index);

                        // Remove the two triangles
                        remove_triangle(tri_index);
                        remove_triangle(neighbour_tri_index);

                        // Add two triangles
                        int first_triangle_id, second_triangle_id;

                        first_triangle_id = add_triangle(pt, neighbour_tri_other_vertex, edge_a_pt);
                        second_triangle_id = add_triangle(pt, neighbour_tri_other_vertex, edge_b_pt);

                        // recursion below
                        flip_bad_edges(first_triangle_id, pt);
                        flip_bad_edges(second_triangle_id, pt);
                    }
                }
            }

            private int get_unique_edge_id()
            {
                int edge_id;
                // get an unique edge id
                if (unique_edgeid_list.Count != 0)
                {
                    edge_id = unique_edgeid_list[0]; // retrive the edge id from the list which stores the id of deleted edges
                    unique_edgeid_list.RemoveAt(0); // remove that id from the unique edge id list
                }
                else
                {
                    edge_id = this.all_edges.Count;
                }
                return edge_id;

            }

            private int get_unique_triangle_id()
            {
                int tri_id;
                // get an unique triangle id
                if (unique_triangleid_list.Count != 0)
                {
                    tri_id = unique_triangleid_list[0]; // retrive the triangle id from the list which stores the id of deleted edges
                    unique_triangleid_list.RemoveAt(0); // remove that id from the unique triangle id list
                }
                else
                {
                    tri_id = this.all_triangles.Count;
                }
                return tri_id;
            }

            private void set_bounding_triangle(List<point_store> all_input_vertices)
            {
                point_store[] x_sorted = all_input_vertices.OrderBy(obj => obj.x).ToArray();
                point_store[] y_sorted = all_input_vertices.OrderBy(obj => obj.y).ToArray();

                // Define bounding triangle
                double max_x, max_y, k;
                max_x = (x_sorted[x_sorted.Length - 1].x - x_sorted[0].x);
                max_y = (y_sorted[y_sorted.Length - 1].y - y_sorted[0].y);
                k = 1000 * Math.Max(max_x, max_y);

                // zeoth _point
                double x_zero, y_zero;
                x_zero = (x_sorted[x_sorted.Length - 1].x + x_sorted[0].x) * 0.5;
                y_zero = (y_sorted[y_sorted.Length - 1].y + y_sorted[0].y) * 0.5;

                // id for the super triangle points
                int pt_count = all_input_vertices.Count;


                // set the vertex
                this.s_p1 = new point_store(pt_count + 1, 0, Math.Round(k / 2.0f), new pslg_datastructure.point2d(pt_count + 1, x_zero, k));
                this.s_p2 = new point_store(pt_count + 2, Math.Round(k / 2.0f), 0.0, new pslg_datastructure.point2d(pt_count + 2, k, y_zero));
                this.s_p3 = new point_store(pt_count + 3, -1 * Math.Round(k / 2.0f), -1 * Math.Round(k / 2.0f), new pslg_datastructure.point2d(pt_count + 3, -k, -k));

                // set the edges
                add_triangle(this.s_p1, this.s_p2, this.s_p3);
            }

            public class point_store
            {
                int _pt_id;
                double _x;
                double _y;
                pslg_datastructure.point2d parent_pt;

                public int pt_id
                {
                    get { return this._pt_id; }
                }

                public double x
                {
                    get { return this._x; }
                }

                public double y
                {
                    get { return this._y; }
                }

                public pslg_datastructure.point2d get_parent_data_type
                {
                    get { return parent_pt; }
                }

                public point_store(int i_pt_id, double i_x, double i_y, pslg_datastructure.point2d pt)
                {
                    // Constructor
                    // set id
                    this._pt_id = i_pt_id;
                    // co-ordinate
                    this._x = i_x;
                    this._y = i_y;
                    // store the output point variable
                    this.parent_pt = pt;
                }

                public bool Equals(point_store other)
                {
                    if (Math.Abs(this._x - other.x) <= eps && Math.Abs(this._y - other.y) <= eps)
                    {
                        return true;
                    }
                    return false;
                }

                // Operators
                public point_store sub(point_store other_pt)
                {
                    double ab_x = this.x - other_pt.x;
                    double ab_y = this.y - other_pt.y;

                    return new point_store(-1, ab_x, ab_y, null);
                }

                public point_store add(point_store other_pt)
                {
                    double ab_x = this.x + other_pt.x;
                    double ab_y = this.y + other_pt.y;

                    return new point_store(-1, ab_x, ab_y, null);
                }

                public double dot(point_store other_pt)
                {
                    return ((this.x * other_pt.x) + (this.y * other_pt.y));
                }

                public double cross(point_store other_pt)
                {
                    return ((this.y * other_pt.x) - (this.x * other_pt.y));
                }

                public point_store mult(double v)
                {
                    return (new point_store(this.pt_id, this.x * v, this.y * v, null));
                }

            }

            public class edge_store
            {
                int _edge_id;
                point_store _start_pt;
                point_store _end_pt;
                triangle_store _left_triangle;
                triangle_store _right_triangle;

                public int edge_id
                {
                    get { return this._edge_id; }
                }

                public point_store start_pt
                {
                    get { return this._start_pt; }
                }

                public point_store end_pt
                {
                    get { return this._end_pt; }
                }

                public edge_store sym
                {
                    get { return new edge_store(this._edge_id, this._end_pt, this._start_pt); }
                }

                public double edge_length
                {
                    get { return Math.Sqrt(Math.Pow(this._start_pt.x - this._end_pt.x, 2) + Math.Pow(this._start_pt.y - this._end_pt.y, 2)); }
                }

                public edge_store(int i_edge_id, point_store i_start_pt, point_store i_end_pt)
                {
                    // Constructor
                    // set id
                    this._edge_id = i_edge_id;
                    // set start and end pt
                    this._start_pt = i_start_pt;
                    this._end_pt = i_end_pt;

                    // set triangles to null
                    this._left_triangle = null;
                    this._right_triangle = null;
                }

                public bool contains_point(point_store the_point)
                {
                    // find whether the point belongs to the triangle
                    if (the_point.Equals(this._start_pt) == true ||
                        the_point.Equals(this._end_pt) == true)
                    {
                        return true;
                    }
                    return false;
                }

                public int other_triangle_id(triangle_store the_triangle)
                {
                    // Function to return the other triangle associated with this edge
                    if (this._left_triangle != null)
                    {
                        if (the_triangle.Equals(this._left_triangle) == true)
                        {
                            if (this._right_triangle == null)
                            {
                                return -1;
                            }
                            return this._right_triangle.tri_id;
                        }
                    }

                    if (this._right_triangle != null)
                    {
                        if (the_triangle.Equals(this._right_triangle) == true)
                        {
                            if (this._left_triangle == null)
                            {
                                return -1;
                            }
                            return this._left_triangle.tri_id;
                        }
                    }

                    return -1;

                }

                public void add_triangle(triangle_store the_triangle)
                {
                    // check whether the input triangle has this edge
                    if (the_triangle.contains_edge(this) == true)
                    {
                        if (rightof(the_triangle.mid_pt, this) == true)
                        {
                            // Add the right triangle
                            this._right_triangle = the_triangle;
                        }
                        else
                        {
                            // Add the left triangle
                            this._left_triangle = the_triangle;
                        }
                    }
                }

                public void remove_triangle(triangle_store the_triangle)
                {
                    // check whether the input triangle has this edge
                    if (the_triangle.contains_edge(this) == true)
                    {
                        if (rightof(the_triangle.mid_pt, this) == true)
                        {
                            // Remove the right triangle
                            this._right_triangle = null;
                        }
                        else
                        {
                            // Remove the left triangle
                            this._left_triangle = null;
                        }
                    }

                }

                private bool ccw(point_store a, point_store b, point_store c)
                {
                    // Computes | a.x a.y  1 |
                    //          | b.x b.y  1 | > 0
                    //          | c.x c.y  1 |
                    return (((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x))) > 0;
                }

                private bool rightof(point_store x, edge_store e)
                {
                    return ccw(x, e.end_pt, e.start_pt);
                }

                public bool Equals(edge_store other)
                {
                    if (other.start_pt.Equals(this._start_pt) == true && other.end_pt.Equals(this._end_pt) == true)
                    {
                        return true;
                    }
                    return false;
                }

                public bool Equals_without_orientation(edge_store other)
                {
                    if ((other.Equals(this) == true) || (other.Equals(this.sym) == true))
                    {
                        return true;
                    }
                    return false;
                }

                public double find_closest_distance(point_store pt)
                {
                    point_store ab = end_pt.sub(start_pt); // ab_x = x2 - x1,  ab_y = y2 - y1  (end_pt - start_pt)
                    double ab_dot = ab.dot(ab); // ab_x*ab_x + ab_y*ab_y

                    // double u=((x3 - x1) * px + (y3 - y1) * py) / ((px*px)+(py*py));
                    double u = (pt.sub(start_pt)).dot(ab) / ab_dot;

                    if (u < 0.0)
                    {
                        u = 0.0;
                    }
                    else if (u > 1.0)
                    {
                        u = 1.0;
                    }

                    // closest point on the edge from the input pt
                    point_store c = start_pt.add(ab.mult(u));
                    // vector connecting the input pt and the pt on line
                    point_store pt_to_c = c.sub(pt);

                    return (pt_to_c.dot(pt_to_c));
                }

                public point_store find_common_pt(edge_store other_edge)
                {
                    // Returns the common point of this edge and the other edges
                    if (this.start_pt.Equals(other_edge.start_pt))
                    {
                        return this.start_pt;
                    }
                    else if (this.start_pt.Equals(other_edge.end_pt))
                    {
                        return this.start_pt;
                    }
                    else if (this.end_pt.Equals(other_edge.start_pt))
                    {
                        return this.end_pt;
                    }
                    else if (this.end_pt.Equals(other_edge.end_pt))
                    {
                        return this.end_pt;
                    }
                    return null;
                }

                public point_store find_other_pt(point_store pt)
                {
                    if (this.start_pt.Equals(pt) == true)
                    {
                        return this.end_pt;
                    }
                    else if (this.end_pt.Equals(pt) == true)
                    {
                        return this.start_pt;
                    }

                    return null;
                }


            }

            public class triangle_store
            {
                int _tri_id;
                point_store _pt1;
                point_store _pt2;
                point_store _pt3;
                point_store _mid_pt;
                edge_store _e1; // from pt1 to pt2
                edge_store _e2; // from pt2 to pt3
                edge_store _e3; // from pt3 to pt1
                point_store _circum_center;
                double _circum_circle_radius;
                double _shortest_edge;

                public int tri_id
                {
                    get { return this._tri_id; }
                }

                public point_store pt1
                {
                    get { return this._pt1; }
                }

                public point_store pt2
                {
                    get { return this._pt2; }
                }

                public point_store pt3
                {
                    get { return this._pt3; }
                }

                public edge_store e1
                {
                    get { return this._e1; }
                }

                public edge_store e2
                {
                    get { return this._e2; }
                }

                public edge_store e3
                {
                    get { return this._e3; }
                }

                public point_store mid_pt
                {
                    // return the mid point
                    get { return this._mid_pt; }
                }

                public point_store circum_center
                {
                    get { return this._circum_center; }
                }

                public double circum_radius
                {
                    get { return this._circum_circle_radius; }
                }

                public double circumradius_shortest_edge_ratio
                {
                    get { return (this._circum_circle_radius / this._shortest_edge); }
                }

                public triangle_store(int i_tri_id, point_store i_p1, point_store i_p2, point_store i_p3, edge_store i_e1, edge_store i_e2, edge_store i_e3)
                {
                    // Constructor
                    // set id
                    this._tri_id = i_tri_id;
                    // set points
                    this._pt1 = i_p1;
                    this._pt2 = i_p2;
                    this._pt3 = i_p3;
                    // set edges
                    this._e1 = i_e1;
                    this._e2 = i_e2;
                    this._e3 = i_e3;
                    // set the mid point
                    this._mid_pt = new point_store(-1, (this._pt1.x + this._pt2.x + this._pt3.x) / 3.0f, (this._pt1.y + this._pt2.y + this._pt3.y) / 3.0f, null);

                    set_incircle();
                    set_shortest_edge();
                }

                private void set_incircle()
                {
                    double dA = (pt1.x * pt1.x) + (pt1.y * pt1.y);
                    double dB = (pt2.x * pt2.x) + (pt2.y * pt2.y);
                    double dC = (pt3.x * pt3.x) + (pt3.y * pt3.y);

                    double aux1 = (dA * (pt3.y - pt2.y) + dB * (pt1.y - pt3.y) + dC * (pt2.y - pt1.y));
                    double aux2 = -(dA * (pt3.x - pt2.x) + dB * (pt1.x - pt3.x) + dC * (pt2.x - pt1.x));
                    double div = (2 * (pt1.x * (pt3.y - pt2.y) + pt2.x * (pt1.y - pt3.y) + pt3.x * (pt2.y - pt1.y)));

                    //Circumcircle
                    double center_x = aux1 / div;
                    double center_y = aux2 / div;

                    this._circum_center = new point_store(-1, center_x, center_y, new pslg_datastructure.point2d(-100, center_x, center_y));
                    this._circum_circle_radius = Math.Sqrt((center_x - pt1.x) * (center_x - pt1.x) + (center_y - pt1.y) * (center_y - pt1.y));
                    //this._ellipse_edge = new point2d(-1, center_x - this._circle_radius, center_y - this._circle_radius);
                    // this._ellipse_edge = new point2d(-1, center_x - 2, center_y - 2);
                }

                private void set_shortest_edge()
                {
                    double edge_len_1 = e1.edge_length;
                    double edge_len_2 = e2.edge_length;
                    double edge_len_3 = e3.edge_length;

                    this._shortest_edge = edge_len_1 < edge_len_2 ? (edge_len_1 < edge_len_3 ? edge_len_1 : edge_len_3) : (edge_len_2 < edge_len_3 ? edge_len_2 : edge_len_3);
                    //x < y ? (x < z ? x : z) : (y < z ? y : z)
                }

                public bool contains_point(point_store the_point)
                {
                    // find whether the point belongs to the triangle
                    if (the_point.Equals(this._pt1) == true ||
                        the_point.Equals(this._pt2) == true ||
                        the_point.Equals(this._pt3) == true)
                    {
                        return true;
                    }
                    return false;
                }

                public point_store get_other_pt(edge_store the_edge)
                {
                    // Returns the third vertex of this triangle which is not part of the given edge
                    if (the_edge.Equals_without_orientation(this.e1) == true)
                    {
                        return (this.e2.find_common_pt(this.e3));
                    }
                    else if (the_edge.Equals_without_orientation(this.e2) == true)
                    {
                        return (this.e3.find_common_pt(this.e1));
                    }
                    else if (the_edge.Equals_without_orientation(this.e3) == true)
                    {
                        return (this.e1.find_common_pt(this.e2));
                    }
                    return null;
                }

                public int get_other_edge_id(point_store pt)
                {
                    if (this.e1.start_pt.Equals(pt) == false && this.e1.end_pt.Equals(pt) == false)
                    {
                        return this.e1.edge_id;
                    }
                    else if (this.e2.start_pt.Equals(pt) == false && this.e2.end_pt.Equals(pt) == false)
                    {
                        return this.e2.edge_id;
                    }
                    else if (this.e3.start_pt.Equals(pt) == false && this.e3.end_pt.Equals(pt) == false)
                    {
                        return this.e3.edge_id;
                    }

                    return -1;
                }


                public bool contains_edge(edge_store the_edge)
                {
                    // find whether the edge belongs to the triangle
                    if (the_edge.Equals_without_orientation(this._e1) == true ||
                        the_edge.Equals_without_orientation(this._e2) == true ||
                        the_edge.Equals_without_orientation(this._e3) == true)
                    {
                        return true;
                    }
                    return false;
                }


                // Tests if a 2D point lies inside this 2D triangle.See Real-Time Collision
                // Detection, chap. 5, p. 206.
                //
                // @param point
                // The point to be tested
                // @return Returns true iff the point lies inside this 2D triangle
                public bool is_point_inside(point_store the_pt)
                {
                    double pab = the_pt.sub(this.pt1).cross(this.pt2.sub(this.pt1));
                    double pbc = the_pt.sub(this.pt2).cross(this.pt3.sub(this.pt2));

                    if (has_same_sign(pab, pbc) == false)
                    {
                        return false;
                    }

                    double pca = the_pt.sub(this.pt3).cross(this.pt1.sub(this.pt3));

                    return has_same_sign(pab, pca);
                }

                private bool has_same_sign(double a, double b)
                {
                    return Math.Sign(a) == Math.Sign(b);
                }


                //// Tests if a given point lies in the circumcircle of this triangle.Let the
                //// triangle ABC appear in counterclockwise (CCW) order. Then when det &gt;
                //// 0, the point lies inside the circumcircle through the three points a, b
                //// and c.If instead det &lt; 0, the point lies outside the circumcircle.
                //// When det = 0, the four points are cocircular.If the triangle is oriented
                //// clockwise (CW) the result is reversed.See Real-Time Collision Detection,
                //// chap. 3, p. 34.

                public bool is_point_in_circumcircle(point_store the_pt)
                {
                    double a11 = pt1.x - the_pt.x;
                    double a21 = pt2.x - the_pt.x;
                    double a31 = pt3.x - the_pt.x;

                    double a12 = pt1.y - the_pt.y;
                    double a22 = pt2.y - the_pt.y;
                    double a32 = pt3.y - the_pt.y;

                    double a13 = (pt1.x - the_pt.x) * (pt1.x - the_pt.x) + (pt1.y - the_pt.y) * (pt1.y - the_pt.y);
                    double a23 = (pt2.x - the_pt.x) * (pt2.x - the_pt.x) + (pt2.y - the_pt.y) * (pt2.y - the_pt.y);
                    double a33 = (pt3.x - the_pt.x) * (pt3.x - the_pt.x) + (pt3.y - the_pt.y) * (pt3.y - the_pt.y);

                    double det = (a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 - a12 * a21 * a33 - a11 * a23 * a32);

                    return ((is_oriented_ccw() == true) ? det > 0.0 : det < 0.0);
                }

                public bool is_point_inside_circumcircle(point_store pt)
                {
                    double d_squared = (pt.x - circum_center.x) * (pt.x - circum_center.x) +
                        (pt.y - circum_center.y) * (pt.y - circum_center.y);
                    return d_squared < (this._circum_circle_radius * this._circum_circle_radius);
                }

                private bool is_oriented_ccw()
                {
                    double a11 = pt1.x - pt3.x;
                    double a21 = pt2.x - pt3.x;

                    double a12 = pt1.y - pt3.y;
                    double a22 = pt2.y - pt3.y;

                    double det = a11 * a22 - a12 * a21;

                    return det > 0.0;
                }


                public bool Equals(triangle_store the_triangle)
                {
                    // find the triangle equals
                    if (this.contains_edge(the_triangle.e1) == true &&
                       this.contains_edge(the_triangle.e2) == true &&
                       this.contains_edge(the_triangle.e3) == true)
                    {
                        return true;
                    }
                    return false;
                }

            }

        }
        #endregion
    }
}
