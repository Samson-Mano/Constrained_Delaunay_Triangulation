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

            // step : 5 Well sized triangle condition
            // parameter 1 => B is the parameter from Chew's first algorithm B = 1, Ruppert's B = Sqrt(2), Chews second algorithm B = Sqrt(5)/2
            double B_var = Math.Sqrt(2);
            // parameter 2 => h is the desired side length of triangle in the triangulation (intiutive user input)
            double h_var = main_mesh.all_triangles.OrderBy(obj => obj.shortest_edge).ToArray()[0].shortest_edge * 3.0f; // change this !!!! to desired size

            // step : 6 Find and queue the bad triangles
            pslg_datastructure.surface_store current_surface = the_surface_data[the_surface_index];
            mesh_store.triangle_store bad_triangle = main_mesh.all_triangles.Find(obj => triangle_angle_size_constraint(current_surface, outter_edges, inner_edges, obj, B_var, h_var) == true);

            while (bad_triangle != null)
            {
                // c_vertex stores the circum_center of bad triangles
                pslg_datastructure.point2d inner_surface_pt = new pslg_datastructure.point2d(-1, bad_triangle.circum_center.x, bad_triangle.circum_center.y);

                // step : 6A Refine the outter edges which are encroched by the failed triangles circum center
                is_encroched = false;
                encroched_segment_single_divide(ref outter_edges, ref inner_surface_pt, h_var, ref is_encroched, ref seed_outter_edge_exists);

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
                    encroched_segment_single_divide(ref inner_edges[j], ref inner_surface_pt, h_var, ref is_encroched, ref seed_inner_edge_exists[j]);

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
                // Find the new bad triangle
                bad_triangle = main_mesh.all_triangles.Find(obj => triangle_angle_size_constraint(current_surface, outter_edges, inner_edges, obj, B_var, h_var) == true);
            }

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
                            // Point where the segement is split
                            pslg_datastructure.point2d split_pt = new pslg_datastructure.point2d(all_nodes.Count + amend_nodes.Count, surface_edges[i].mid_pt.x, surface_edges[i].mid_pt.y);

                            // Check for acute angle connected segement
                            // This sequence is to avoid non-convergence in spliting small angled segement infinite times
                            // Find the two segments connected to the point
                            List<pslg_datastructure.edge2d> two_segments = surface_edges.FindAll(obj => obj.start_pt.Equals(pt) || obj.end_pt.Equals(pt));

                            if (two_segments.Count != 0)
                            {
                                // if (two_segments.Count >2) // Some thing is terribly wrong
                                pslg_datastructure.point2d other_pt1 = two_segments[0].start_pt.Equals(pt) == true ? two_segments[0].end_pt : two_segments[0].start_pt;
                                pslg_datastructure.point2d other_pt2 = two_segments[two_segments.Count - 1].end_pt.Equals(pt) == true ? two_segments[two_segments.Count - 1].start_pt : two_segments[two_segments.Count - 1].end_pt;

                                if (surface_edges[i].vertex_exists(other_pt1) == true || surface_edges[i].vertex_exists(other_pt2) == true)
                                {

                                    pslg_datastructure.point2d apex_pt, non_apex_vertex_seg1, non_apex_vertex_seg2;
                                    // The edge and the point encroched shares the same vertex (call it apex vertex)
                                    double split_length, t_param;
                                    // I dont really understand why Ruppert & JR Shewchunk using concentric circles of radii of two squares
                                    // By using a split lenght of smallest edge we end up with a isoceles triangle
                                    if (surface_edges[i].vertex_exists(other_pt1) == true)
                                    {
                                        apex_pt = other_pt1;
                                        // Find the other vertex of the two segments
                                        non_apex_vertex_seg1 = surface_edges[i].start_pt.Equals(apex_pt) == true ? surface_edges[i].end_pt : surface_edges[i].start_pt;
                                        non_apex_vertex_seg2 = two_segments[0].start_pt.Equals(apex_pt) == true ? two_segments[0].end_pt : two_segments[0].start_pt;
                                        // edge length of the smallest segment is chosen
                                        split_length = two_segments[0].edge_length;
                                    }
                                    else
                                    {
                                        apex_pt = other_pt2;
                                        // Find the other vertex of the two segments
                                        non_apex_vertex_seg1 = surface_edges[i].start_pt.Equals(apex_pt) == true ? surface_edges[i].end_pt : surface_edges[i].start_pt;
                                        non_apex_vertex_seg2 = two_segments[0].start_pt.Equals(apex_pt) == true ? two_segments[two_segments.Count - 1].end_pt : two_segments[two_segments.Count - 1].start_pt;

                                        // edge length of the smallest segment is chosen
                                        split_length = two_segments[two_segments.Count - 1].edge_length;

                                    }

                                    // Check whether the angle made is less than 30 degreee
                                    if (Form1.the_static_class.GetAngle(non_apex_vertex_seg1.x, non_apex_vertex_seg1.y, apex_pt.x, apex_pt.y, non_apex_vertex_seg2.x, non_apex_vertex_seg2.y) < 0.5236)
                                    {

                                        t_param = split_length / surface_edges[i].edge_length;
                                        if (t_param > 0.5f)
                                        {
                                            double split_x = apex_pt.x * (1 - t_param) + (non_apex_vertex_seg1.x * t_param);
                                            double split_y = apex_pt.y * (1 - t_param) + (non_apex_vertex_seg1.y * t_param);

                                            split_pt = new pslg_datastructure.point2d(all_nodes.Count + amend_nodes.Count, split_x, split_y);
                                        }
                                    }

                                }
                            }

                            // create two segements from one segment by splitting at the mid point
                            pslg_datastructure.edge2d seg_1 = new pslg_datastructure.edge2d(surface_edges[i].edge_id, surface_edges[i].start_pt, split_pt); // create a segment 1 with start_pt to mid_pt
                            pslg_datastructure.edge2d seg_2 = new pslg_datastructure.edge2d(surface_edges.Count, split_pt, surface_edges[i].end_pt); // create a segment 2 with mid_pt to end_pt
                            amend_nodes.Add(new pslg_datastructure.point2d(all_nodes.Count + amend_nodes.Count, split_pt.x, split_pt.y));

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

        public static void encroched_segment_single_divide(ref List<pslg_datastructure.edge2d> surface_edges, ref pslg_datastructure.point2d circumcenter_pt, double meansize, ref bool to_continue, ref bool seed_exists)
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
                        // Find whether the input vertex is attached to any of the edges
                        pslg_datastructure.point2d c_pt = circumcenter_pt;
                        if (surface_edges.Exists(obj => obj.start_pt.Equals(c_pt) || obj.end_pt.Equals(c_pt)) == true)
                        {
                            // If the input vertex attached to any of the edges then find whether the encroched edge and vertex edge forms a small angle
                            // Find the two segments connected to the point
                            List<pslg_datastructure.edge2d> two_segments = surface_edges.FindAll(obj => obj.start_pt.Equals(c_pt) || obj.end_pt.Equals(c_pt));

                            // if (two_segments.Count >2) // Some thing is terribly wrong
                            pslg_datastructure.point2d other_pt1 = two_segments[0].start_pt.Equals(c_pt) == true ? two_segments[0].end_pt : two_segments[0].start_pt;
                            pslg_datastructure.point2d other_pt2 = two_segments[two_segments.Count - 1].end_pt.Equals(c_pt) == true ? two_segments[two_segments.Count - 1].start_pt : two_segments[two_segments.Count - 1].end_pt;


                            if (surface_edges[i].vertex_exists(other_pt1) == true || surface_edges[i].vertex_exists(other_pt2) == true)
                            {
                                // Find the apex point
                                pslg_datastructure.point2d apex_pt, non_apex_vertex_seg1, non_apex_vertex_seg2;
                                // The edge and the point encroched shares the same vertex (call it apex vertex)
                                double split_length, t_param;
                                // I dont really understand why Ruppert & JR Shewchunk using concentric circles of radii of two squares
                                // By using a split lenght of smallest edge we end up with a isoceles triangle
                                if (surface_edges[i].vertex_exists(other_pt1) == true)
                                {
                                    apex_pt = other_pt1;
                                    // Find the other vertex of the two segments
                                    non_apex_vertex_seg1 = surface_edges[i].start_pt.Equals(apex_pt) == true ? surface_edges[i].end_pt : surface_edges[i].start_pt;
                                    non_apex_vertex_seg2 = two_segments[0].start_pt.Equals(apex_pt) == true ? two_segments[0].end_pt : two_segments[0].start_pt;
                                    // edge length of the smallest segment is chosen
                                    split_length = two_segments[0].edge_length;
                                }
                                else
                                {
                                    apex_pt = other_pt2;
                                    // Find the other vertex of the two segments
                                    non_apex_vertex_seg1 = surface_edges[i].start_pt.Equals(apex_pt) == true ? surface_edges[i].end_pt : surface_edges[i].start_pt;
                                    non_apex_vertex_seg2 = two_segments[0].start_pt.Equals(apex_pt) == true ? two_segments[two_segments.Count - 1].end_pt : two_segments[two_segments.Count - 1].start_pt;

                                    // edge length of the smallest segment is chosen
                                    split_length = two_segments[two_segments.Count - 1].edge_length;

                                }

                                // Check whether the angle made is less than 30 degreee
                                if (Form1.the_static_class.GetAngle(non_apex_vertex_seg1.x, non_apex_vertex_seg1.y, apex_pt.x, apex_pt.y, non_apex_vertex_seg2.x, non_apex_vertex_seg2.y) < 0.5236)
                                {
                                    // Yes the angle is less than 30 degree making it a small angle case
                                    t_param = split_length / surface_edges[i].edge_length;
                                    if (t_param > 0.5f)
                                    {
                                        double split_x = apex_pt.x * (1 - t_param) + (non_apex_vertex_seg1.x * t_param);
                                        double split_y = apex_pt.y * (1 - t_param) + (non_apex_vertex_seg1.y * t_param);

                                        circumcenter_pt = new pslg_datastructure.point2d(-1, split_x, split_y);

                                        pslg_datastructure.edge2d tseg_1 = new pslg_datastructure.edge2d(surface_edges[i].edge_id, surface_edges[i].start_pt, circumcenter_pt); // create a segment 1 with start_pt to mid_pt
                                        pslg_datastructure.edge2d tseg_2 = new pslg_datastructure.edge2d(surface_edges.Count, circumcenter_pt, surface_edges[i].end_pt); // create a segment 2 with mid_pt to end_pt

                                        surface_edges.RemoveAt(i); // remove the edge which encroches

                                        // Insert the newly created two element edge list to the main list (at the removed index)
                                        surface_edges.Insert(i, tseg_1);
                                        surface_edges.Insert(i + 1, tseg_2);
                                        to_continue = true;

                                        // break is very important to avoid no longer continuing the search for the pt encroching i_th index segment
                                        break;
                                    }
                                    else
                                    {
                                        // the ratio of length is not more than half which means the longest edge is twice the length of smallest edge
                                        goto half_split;
                                    }
                                }
                                else
                                {
                                    // Angle is not less than 30 degree so need to impose the split
                                    goto half_split;
                                }
                            }
                            else
                            {
                                // No apex point
                                // which means no relation to the circum center point so bisect this edge into half
                                // create two segements from one segment by splitting at the mid point
                                goto half_split;

                            }
                        }

                    half_split:;
                        // the vertex point is actualy a circle center
                        // create two segements from one segment by splitting at the mid point
                        pslg_datastructure.edge2d seg_1 = new pslg_datastructure.edge2d(surface_edges[i].edge_id, surface_edges[i].start_pt, surface_edges[i].mid_pt); // create a segment 1 with start_pt to mid_pt
                        pslg_datastructure.edge2d seg_2 = new pslg_datastructure.edge2d(surface_edges.Count, surface_edges[i].mid_pt, surface_edges[i].end_pt); // create a segment 2 with mid_pt to end_pt
                        circumcenter_pt = new pslg_datastructure.point2d(-1, surface_edges[i].mid_pt.x, surface_edges[i].mid_pt.y);

                        surface_edges.RemoveAt(i); // remove the edge which encroches

                        // Insert the newly created two element edge list to the main list (at the removed index)
                        surface_edges.Insert(i, seg_1);
                        surface_edges.Insert(i + 1, seg_2);
                        to_continue = true;

                        // break is very important to avoid no longer continuing the search for the pt encroching i_th index segment
                        break;
                    }
                }
            }
        }

        public static bool triangle_angle_size_constraint(pslg_datastructure.surface_store the_surface, List<pslg_datastructure.edge2d> outter_edges, List<pslg_datastructure.edge2d>[] inner_edges, mesh_store.triangle_store graded_triangle, double B, double h)
        {
            if (the_surface.pointinsurface(graded_triangle.shrunk_vertices[0].x, graded_triangle.shrunk_vertices[0].y) == true &&
                        the_surface.pointinsurface(graded_triangle.shrunk_vertices[1].x, graded_triangle.shrunk_vertices[1].y) == true &&
                        the_surface.pointinsurface(graded_triangle.shrunk_vertices[2].x, graded_triangle.shrunk_vertices[2].y) == true)
            {
                if (graded_triangle.circumradius_shortest_edge_ratio > B)
                {
                    // condition 1: B parameter => A triangle is well-shaped if all its angles are greater than or equal to 30 degrees
                    /* ---------------- Small Angle Input Case ----------------------------------- Very Expensive
                    // Test whether the small angle is due to input from user
                    pslg_datastructure.edge2d temp_e1 = new pslg_datastructure.edge2d(-1, graded_triangle.e1.start_pt.get_parent_data_type, graded_triangle.e1.end_pt.get_parent_data_type);
                    pslg_datastructure.edge2d temp_e2 = new pslg_datastructure.edge2d(-1, graded_triangle.e2.start_pt.get_parent_data_type, graded_triangle.e2.end_pt.get_parent_data_type);
                    pslg_datastructure.edge2d temp_e3 = new pslg_datastructure.edge2d(-1, graded_triangle.e3.start_pt.get_parent_data_type, graded_triangle.e3.end_pt.get_parent_data_type);

                    int e1_index, e2_index, e3_index;
                    e1_index = outter_edges.FindIndex(obj => obj.Equals_without_orientation(temp_e1));
                    e2_index = outter_edges.FindIndex(obj => obj.Equals_without_orientation(temp_e2));
                    e3_index = outter_edges.FindIndex(obj => obj.Equals_without_orientation(temp_e3));

                    if (e1_index != -1 || e2_index != -1 || e3_index != -1)
                    {
                        if (e1_index != -1 && e2_index != -1)
                        {
                            return false;
                        }
                        else if (e2_index != -1 && e3_index != -1)
                        {
                            return false;
                        }
                        else if (e3_index != -1 && e1_index != -1)
                        {
                            return false;
                        }
                    }

                    for (int i = 0; i < inner_edges.Length - 1; i++)
                    {
                        e1_index = inner_edges[i].FindIndex(obj => obj.Equals_without_orientation(temp_e1));
                        e2_index = inner_edges[i].FindIndex(obj => obj.Equals_without_orientation(temp_e2));
                        e3_index = inner_edges[i].FindIndex(obj => obj.Equals_without_orientation(temp_e3));

                        if (e1_index != -1 || e2_index != -1 || e3_index != -1)
                        {
                            if (e1_index != -1 && e2_index != -1)
                            {
                                return false;
                            }
                            else if (e2_index != -1 && e3_index != -1)
                            {
                                return false;
                            }
                            else if (e3_index != -1 && e1_index != -1)
                            {
                                return false;
                            }
                        }
                    }
                    */
                    return true;
                }
                else if (graded_triangle.circum_radius > h && graded_triangle.shortest_edge > h)
                {
                    // condition 2: h parameter => A triangle is well-sized if it satisfies a user - supplied grading function
                    return true;
                }
                else
                {
                    // Conditions are not met
                    return false;
                }
            }
            else
            {
                return false;
            }

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
                point_store temp_pt = new point_store(this.all_points.Count, parent_onemore_point.x, parent_onemore_point.y, parent_onemore_point);
                // Add the point to local list
                this._all_points.Add(temp_pt);
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
                    // collect the edges of the triangle
                    edge_store tri_edge_a = this.all_triangles[containing_triangle_index].e1;
                    edge_store tri_edge_b = this.all_triangles[containing_triangle_index].e2;
                    edge_store tri_edge_c = this.all_triangles[containing_triangle_index].e3;

                    // remove the single triangle
                    remove_triangle(containing_triangle_index);

                    // add the three triangles
                    int[] triangle_id = new int[3];
                    triangle_id = add_three_triangles(pt, tri_edge_a, tri_edge_b, tri_edge_c);

                    // Flip the bad triangles recursively
                    flip_bad_edges(triangle_id[0], pt);
                    flip_bad_edges(triangle_id[1], pt);
                    flip_bad_edges(triangle_id[2], pt);
                }
                else
                {
                    // Point lies on the edge
                    // Find the edge which is closest to the pt
                    int the_edge_id = this.all_edges.Find(obj => obj.test_point_on_line(pt)).edge_id;

                    int first_tri_index = this.all_edges[the_edge_id].left_triangle.tri_id;
                    int second_tri_index = this.all_edges[the_edge_id].right_triangle.tri_id;

                    // collect the other two edges of first triangle
                    edge_store[] first_tri_other_two_edge = new edge_store[2];
                    first_tri_other_two_edge = this.all_triangles[first_tri_index].get_other_two_edges(this.all_edges[the_edge_id]);

                    // collect the other two edges of second triangle
                    edge_store[] second_tri_other_two_edge = new edge_store[2];
                    second_tri_other_two_edge = this.all_triangles[second_tri_index].get_other_two_edges(this.all_edges[the_edge_id]);

                    // Remove the common edge
                    unique_edgeid_list.Add(this.all_edges[the_edge_id].edge_id);
                    // this._all_edges.RemoveAt(the_edge_index);
                    this._all_edges[the_edge_id].remove_edge();

                    // remove the two triangle
                    remove_triangle(first_tri_index);
                    remove_triangle(second_tri_index);

                    // add the three triangles
                    int[] triangle_id = new int[4];
                    triangle_id = add_four_triangles(pt, first_tri_other_two_edge[0], first_tri_other_two_edge[1], second_tri_other_two_edge[0], second_tri_other_two_edge[1]);

                    // Flip the bad triangles recursively
                    flip_bad_edges(triangle_id[0], pt);
                    flip_bad_edges(triangle_id[1], pt);
                    flip_bad_edges(triangle_id[2], pt);
                    flip_bad_edges(triangle_id[3], pt);
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
                  

                    if (the_surface.pointinsurface(the_tri.shrunk_vertices[0].x, the_tri.shrunk_vertices[0].y) == false ||
                        the_surface.pointinsurface(the_tri.shrunk_vertices[1].x, the_tri.shrunk_vertices[1].y) == false ||
                        the_surface.pointinsurface(the_tri.shrunk_vertices[2].x, the_tri.shrunk_vertices[2].y) == false)
                    {
                        // continue because the face is not in surface
                        continue;
                    }

                    pslg_datastructure.triangle2d temp_tri = new pslg_datastructure.triangle2d(i, the_tri.pt1.get_parent_data_type, the_tri.pt2.get_parent_data_type, the_tri.pt3.get_parent_data_type);
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
                int edge_index1 = this.all_triangles[tri_index].e1.edge_id;
                int edge_index2 = this.all_triangles[tri_index].e2.edge_id;
                int edge_index3 = this.all_triangles[tri_index].e3.edge_id;

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
                this._all_triangles[tri_index].remove_triangle();
            }

            public int[] add_three_triangles(point_store new_pt, edge_store edge_a, edge_store edge_b, edge_store edge_c)
            {
                // Add three new edges from the new point
                int[] edge_indices = new int[3];
                // Add Edge 1
                edge_indices[0] = add_edge(new_pt, edge_a.start_pt);

                // Add Edge 2
                edge_indices[1] = add_edge(new_pt, edge_b.start_pt);

                // Add Edge 3
                edge_indices[2] = add_edge(new_pt, edge_c.start_pt);
                //_________________________________________________________________________________

                // Add three triangles
                int[] output_indices = new int[3];
                // Add First triangle
                output_indices[0] = add_triangle(new_pt, edge_a.start_pt, edge_a.end_pt, this.all_edges[edge_indices[0]], edge_a, this.all_edges[edge_indices[1]].sym);

                // Add Second triangle
                output_indices[1] = add_triangle(new_pt, edge_b.start_pt, edge_b.end_pt, this.all_edges[edge_indices[1]], edge_b, this.all_edges[edge_indices[2]].sym);

                // Add Third triangle
                output_indices[2] = add_triangle(new_pt, edge_c.start_pt, edge_c.end_pt, this.all_edges[edge_indices[2]], edge_c, this.all_edges[edge_indices[0]].sym);
                //_________________________________________________________________________________

                // Add the triangle details to the edge
                // Edge 1
                this._all_edges[edge_indices[0]].add_triangle(this.all_triangles[output_indices[0]]);
                this._all_edges[edge_indices[0]].add_triangle(this.all_triangles[output_indices[2]]);
                // Edge 2
                this._all_edges[edge_indices[1]].add_triangle(this.all_triangles[output_indices[0]]);
                this._all_edges[edge_indices[1]].add_triangle(this.all_triangles[output_indices[1]]);
                // Edge 3
                this._all_edges[edge_indices[2]].add_triangle(this.all_triangles[output_indices[1]]);
                this._all_edges[edge_indices[2]].add_triangle(this.all_triangles[output_indices[2]]);

                // Edge a
                this._all_edges[edge_a.edge_id].add_triangle(this.all_triangles[output_indices[0]]);
                // Edge b
                this._all_edges[edge_b.edge_id].add_triangle(this.all_triangles[output_indices[1]]);
                // Edge c
                this._all_edges[edge_c.edge_id].add_triangle(this.all_triangles[output_indices[2]]);

                //_________________________________________________________________________________

                return output_indices;
            }

            public int[] add_four_triangles(point_store new_pt, edge_store edge_a, edge_store edge_b, edge_store edge_c, edge_store edge_d)
            {
                // Add four new edges from the new point
                int[] edge_indices = new int[4];
                // Add Edge 1
                edge_indices[0] = add_edge(new_pt, edge_a.start_pt);

                // Add Edge 2
                edge_indices[1] = add_edge(new_pt, edge_b.start_pt);

                // Add Edge 3
                edge_indices[2] = add_edge(new_pt, edge_c.start_pt);

                // Add Edge 4
                edge_indices[3] = add_edge(new_pt, edge_d.start_pt);
                //_________________________________________________________________________________

                // Add four triangles
                int[] output_indices = new int[4];
                // Add First triangle
                output_indices[0] = add_triangle(new_pt, edge_a.start_pt, edge_a.end_pt, this.all_edges[edge_indices[0]], edge_a, this.all_edges[edge_indices[1]].sym);

                // Add Second triangle
                output_indices[1] = add_triangle(new_pt, edge_b.start_pt, edge_b.end_pt, this.all_edges[edge_indices[1]], edge_b, this.all_edges[edge_indices[2]].sym);

                // Add Third triangle
                output_indices[2] = add_triangle(new_pt, edge_c.start_pt, edge_c.end_pt, this.all_edges[edge_indices[2]], edge_c, this.all_edges[edge_indices[3]].sym);

                // Add Fourth triangle
                output_indices[3] = add_triangle(new_pt, edge_d.start_pt, edge_d.end_pt, this.all_edges[edge_indices[3]], edge_d, this.all_edges[edge_indices[0]].sym);
                //_________________________________________________________________________________

                // Add the triangle details to the edge
                // Edge 1
                this._all_edges[edge_indices[0]].add_triangle(this.all_triangles[output_indices[0]]);
                this._all_edges[edge_indices[0]].add_triangle(this.all_triangles[output_indices[3]]);
                // Edge 2
                this._all_edges[edge_indices[1]].add_triangle(this.all_triangles[output_indices[0]]);
                this._all_edges[edge_indices[1]].add_triangle(this.all_triangles[output_indices[1]]);
                // Edge 3
                this._all_edges[edge_indices[2]].add_triangle(this.all_triangles[output_indices[1]]);
                this._all_edges[edge_indices[2]].add_triangle(this.all_triangles[output_indices[2]]);
                // Edge 4
                this._all_edges[edge_indices[3]].add_triangle(this.all_triangles[output_indices[2]]);
                this._all_edges[edge_indices[3]].add_triangle(this.all_triangles[output_indices[3]]);

                // Edge a
                this._all_edges[edge_a.edge_id].add_triangle(this.all_triangles[output_indices[0]]);
                // Edge b
                this._all_edges[edge_b.edge_id].add_triangle(this.all_triangles[output_indices[1]]);
                // Edge c
                this._all_edges[edge_c.edge_id].add_triangle(this.all_triangles[output_indices[2]]);
                // Edge d
                this._all_edges[edge_d.edge_id].add_triangle(this.all_triangles[output_indices[3]]);
                //_________________________________________________________________________________

                return output_indices;
            }

            public int[] add_two_triangles(point_store new_pt, edge_store tri_a_edge_0, edge_store tri_a_edge_1, edge_store tri_b_edge_0, edge_store tri_b_edge_1)
            {
                // add the only new common edge
                int edge_index;
                // Add Edge 
                edge_index = add_edge(tri_a_edge_1.start_pt, tri_b_edge_1.start_pt);
                //_________________________________________________________________________________

                // Add two triangles
                int[] output_indices = new int[2];
                // Add First triangle
                output_indices[0] = add_triangle(tri_a_edge_1.start_pt, tri_b_edge_0.start_pt, this.all_edges[edge_index].sym.start_pt, tri_a_edge_1, tri_b_edge_0, this.all_edges[edge_index].sym);

                // Add Second triangle
                output_indices[1] = add_triangle(tri_b_edge_1.start_pt, tri_a_edge_0.start_pt, this.all_edges[edge_index].start_pt, tri_b_edge_1, tri_a_edge_0, this.all_edges[edge_index]);
                //_________________________________________________________________________________

                // Add the triangle details to the edge
                // Common Edge
                this._all_edges[edge_index].add_triangle(this.all_triangles[output_indices[0]]);
                this._all_edges[edge_index].add_triangle(this.all_triangles[output_indices[1]]);

                // Edge a
                this._all_edges[tri_a_edge_1.edge_id].add_triangle(this.all_triangles[output_indices[0]]);
                this._all_edges[tri_b_edge_0.edge_id].add_triangle(this.all_triangles[output_indices[0]]);
                // Edge b
                this._all_edges[tri_b_edge_1.edge_id].add_triangle(this.all_triangles[output_indices[1]]);
                this._all_edges[tri_a_edge_0.edge_id].add_triangle(this.all_triangles[output_indices[1]]);
                //_________________________________________________________________________________

                return output_indices;
            }

            private int add_edge(point_store pt1, point_store pt2)
            {
                // Add Edge
                int edge_index = get_unique_edge_id();
                if (edge_index > this.all_edges.Count - 1)
                {
                    // new edge is added
                    edge_store temp_edge = new edge_store();
                    temp_edge.add_edge(edge_index, pt1, pt2);
                    this._all_edges.Add(temp_edge);
                }
                else
                {
                    // existing edge is revised (previously deleted edge is now filled)
                    this._all_edges[edge_index].add_edge(edge_index, pt1, pt2);
                }
                return edge_index;
            }

            private int add_triangle(point_store pt1, point_store pt2, point_store pt3, edge_store e1, edge_store e2, edge_store e3)
            {
                // Add Triangle
                int tri_index = get_unique_triangle_id();
                if (tri_index > this.all_triangles.Count - 1)
                {
                    // new triangle is added
                    triangle_store temp_tri = new triangle_store();
                    temp_tri.add_triangle(tri_index, pt1, pt2, pt3, e1, e2, e3);
                    this._all_triangles.Add(temp_tri);
                }
                else
                {
                    // existing triangle is revised (previously deleted triangle is now filled)
                    this._all_triangles[tri_index].add_triangle(tri_index, pt1, pt2, pt3, e1, e2, e3);
                }
                return tri_index;
            }

            private void flip_bad_edges(int tri_index, point_store pt)
            {
                //flip recursively for the new pair of triangles
                //
                //           pl                    pl
                //          /||\                  /  \
                //       al/ || \bl            al/    \bl
                //        /  ||  \              /      \
                //       /  a||b  \    flip    /___a____\
                //     p0\   ||   /pt   =>   p0\---b----/pt
                //        \  ||  /              \      /
                //       ar\ || /br            ar\    /br
                //          \||/                  \  /
                //           p2                    p2
                //
                // find the edge of this triangle whihc does not contain pt
                int common_edge_index = this.all_triangles[tri_index].get_other_edge_id(pt);

                // Find the neighbour tri index
                int neighbour_tri_index = this.all_edges[common_edge_index].other_triangle_id(this.all_triangles[tri_index]);

                // legalize only if the triangle has a neighbour
                if (neighbour_tri_index != -1)
                {
                    // check whether the newly added pt is inside the neighbour triangle
                    if (this.all_triangles[neighbour_tri_index].is_point_in_circumcircle(pt) == true)
                    {
                        // collect the other two edges of the triangle
                        edge_store[] tri_other_two_edges = new edge_store[2];
                        tri_other_two_edges = this.all_triangles[tri_index].get_other_two_edges(this.all_edges[common_edge_index]);

                        // collect the other two edges of the neighbour triangle
                        edge_store[] neighbour_tri_other_two_edges = new edge_store[2];
                        neighbour_tri_other_two_edges = this.all_triangles[neighbour_tri_index].get_other_two_edges(this.all_edges[common_edge_index]);

                        // Remove the common edge
                        unique_edgeid_list.Add(this.all_edges[common_edge_index].edge_id);
                        // this._all_edges.RemoveAt(common_edge_index);
                        this._all_edges[common_edge_index].remove_edge();

                        // Remove the two triangles
                        remove_triangle(tri_index);
                        remove_triangle(neighbour_tri_index);

                        // Add new two triangles
                        int[] triangle_id = new int[2];
                        triangle_id = add_two_triangles(pt, tri_other_two_edges[0], tri_other_two_edges[1], neighbour_tri_other_two_edges[0], neighbour_tri_other_two_edges[1]);

                        // recursion below
                        flip_bad_edges(triangle_id[0], pt);
                        flip_bad_edges(triangle_id[1], pt);
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
                this.s_p1 = new point_store(pt_count + 1, 0, Math.Round(k / 2.0f), new pslg_datastructure.point2d(pt_count + 1, x_zero, -k));
                this.s_p2 = new point_store(pt_count + 2, Math.Round(k / 2.0f), 0.0, new pslg_datastructure.point2d(pt_count + 2, k, y_zero));
                this.s_p3 = new point_store(pt_count + 3, -1 * Math.Round(k / 2.0f), -1 * Math.Round(k / 2.0f), new pslg_datastructure.point2d(pt_count + 3, -k, k));

                // Add three new edges for the triangle from the three new point
                int[] edge_indices = new int[3];
                // Add Edge 1
                edge_indices[0] = add_edge(this.s_p1, this.s_p2);

                // Add Edge 2
                edge_indices[1] = add_edge(this.s_p2, this.s_p3);

                // Add Edge 3
                edge_indices[2] = add_edge(this.s_p3, this.s_p1);
                //_________________________________________________________________________________

                // Create the super triangle
                int sup_tri_index;
                sup_tri_index = add_triangle(this.s_p1, this.s_p2, this.s_p3, this.all_edges[edge_indices[0]], this.all_edges[edge_indices[1]], this.all_edges[edge_indices[2]]);

                // Add the triangle details to the outter edges
                // Edge 1
                this.all_edges[edge_indices[0]].add_triangle(this.all_triangles[sup_tri_index]);
                // Edge 2
                this.all_edges[edge_indices[1]].add_triangle(this.all_triangles[sup_tri_index]);
                // Edge 3
                this.all_edges[edge_indices[2]].add_triangle(this.all_triangles[sup_tri_index]);

                // Super triangle creation complete
                //_________________________________________________________________________________
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
                edge_store _sym;

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
                    get { return this._sym; }
                }

                public triangle_store left_triangle
                {
                    get { return this._left_triangle; }
                }

                public triangle_store right_triangle
                {
                    get { return this._right_triangle; }
                }

                public double edge_length
                {
                    get { return Math.Sqrt(Math.Pow(this._start_pt.x - this._end_pt.x, 2) + Math.Pow(this._start_pt.y - this._end_pt.y, 2)); }
                }

                public edge_store()
                {
                    // Empty Constructor
                }

                public void add_edge(int i_edge_id, point_store i_start_pt, point_store i_end_pt)
                {
                    // set id
                    this._edge_id = i_edge_id;
                    // set start and end pt
                    this._start_pt = i_start_pt;
                    this._end_pt = i_end_pt;

                    // set triangles to null
                    this._left_triangle = null;
                    this._right_triangle = null;
                    //______________________________________________________________

                    // Set the symmetrical edge
                    // set id
                    this._sym = new edge_store();
                    this._sym._edge_id = i_edge_id;
                    // set end and start pt
                    this._sym._start_pt = i_end_pt;
                    this._sym._end_pt = i_start_pt;

                    // set triangles to null
                    this._sym._left_triangle = null;
                    this._sym._right_triangle = null;
                    this._sym._sym = this;
                }

                public void remove_edge()
                {
                    // set id
                    this._edge_id = -1;
                    // remove start and end pt
                    this._start_pt = null;
                    this._end_pt = null;

                    // remove triangles to null
                    this._left_triangle = null;
                    this._right_triangle = null;
                    //______________________________________________________________

                    // Remove the symmetrical edge
                    // set id
                    this._sym._edge_id = -1;
                    // remove end and start pt
                    this._sym._start_pt = null;
                    this._sym._end_pt = null;

                    // remove triangles to null
                    this._sym._left_triangle = null;
                    this._sym._right_triangle = null;
                }

                public int other_triangle_id(triangle_store the_triangle)
                {
                    if (the_triangle.tri_id != -1)
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
                    }
                    return -1;
                }

                public void add_triangle(triangle_store the_triangle)
                {
                    // check whether the input triangle has this edge
                    //if (the_triangle.contains_edge(this) == true)
                    //{
                    if (rightof(the_triangle.mid_pt, this) == true)
                    {
                        // Add the right triangle
                        this._right_triangle = the_triangle;
                        this._sym._left_triangle = the_triangle;
                    }
                    else
                    {
                        // Add the left triangle
                        this._left_triangle = the_triangle;
                        this._sym._right_triangle = the_triangle;
                    }
                    //}
                }

                public void remove_triangle(triangle_store the_triangle)
                {
                    // check whether the input triangle has this edge
                    //if (the_triangle.contains_edge(this) == true)
                    //{
                    if (rightof(the_triangle.mid_pt, this) == true)
                    {
                        // Remove the right triangle
                        this._right_triangle = null;
                        this._sym._left_triangle = null;
                    }
                    else
                    {
                        // Remove the left triangle
                        this._left_triangle = null;
                        this._sym._right_triangle = null;
                    }
                    //}
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
                    if (this.edge_id != -1)
                    {
                        if ((other.Equals(this) == true) || (other.Equals(this.sym) == true))
                        {
                            return true;
                        }
                    }
                    return false;
                }


                public bool test_point_on_line(point_store pt)
                {
                    bool rslt = false;
                    // Step: 1 Find the cross product
                    double dxc, dyc; // Vector 1 Between given point and first point of the line
                    dxc = pt.x - start_pt.x;
                    dyc = pt.y - start_pt.y;

                    double Threshold = edge_length * 0.01;

                    double dx1, dy1; // Vector 2 Between the second and first point of the line
                    dx1 = end_pt.x - start_pt.x;
                    dy1 = end_pt.y - start_pt.y;

                    double crossprd;
                    crossprd = (dxc * dy1) - (dyc * dx1); // Vector cross product

                    if (Math.Abs(crossprd) <= Threshold) // Check whether the cross product is within the threshold (other wise Not on the line)
                    {
                        if (Math.Abs(dx1) >= Math.Abs(dy1)) // The line is more horizontal or <= 45 degrees
                        {
                            rslt = (dx1 > 0 ? (start_pt.x < pt.x && pt.x < end_pt.x ? true : false) :
                                              (end_pt.x < pt.x && pt.x < start_pt.x ? true : false));
                        }
                        else // line is more vertical
                        {
                            rslt = (dy1 > 0 ? (start_pt.y < pt.y && pt.y < end_pt.y ? true : false) :
                                                (end_pt.y < pt.y && pt.y < start_pt.y ? true : false));
                        }
                    }
                    return rslt;
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
                public point_store[] shrunk_vertices { get; } = new point_store[3];

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

                public double circumradius_shortest_edge_ratio
                {
                    get { return (this._circum_circle_radius / this._shortest_edge); }
                }

                public double circum_radius
                {
                    get { return this._circum_circle_radius; }
                }

                public double shortest_edge
                {
                    get { return this._shortest_edge; }
                }

                public triangle_store()
                {
                    // Empty Constructor
                }

                public void add_triangle(int i_tri_id, point_store i_p1, point_store i_p2, point_store i_p3, edge_store i_e1, edge_store i_e2, edge_store i_e3)
                {
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
                    set_shrunk_vertices();
                }

                public void remove_triangle()
                {
                    // set id
                    this._tri_id = -1;
                    // set points
                    this._pt1 = null;
                    this._pt2 = null;
                    this._pt3 = null;
                    // set edges
                    this._e1 = null;
                    this._e2 = null;
                    this._e3 = null;
                    // set the mid point
                    this._mid_pt = null;

                    // in center
                    this._circum_center = null;
                    this._circum_circle_radius = 0.0f;

                    // shortest edge
                    this._shortest_edge = 0.0f;

                    //remove shrunk vertices
                    this.shrunk_vertices[0] = null;
                    this.shrunk_vertices[1] = null;
                    this.shrunk_vertices[2] = null;
                }

                private void set_shrunk_vertices()
                {
                    double shrink_factor = 0.98;
                    point_store pt1 = new point_store(-1, this._mid_pt.x * (1 - shrink_factor) + (this.pt1.x * shrink_factor), this._mid_pt.y * (1 - shrink_factor) + (this.pt1.y * shrink_factor), null);
                    point_store pt2 = new point_store(-1, this._mid_pt.x * (1 - shrink_factor) + (this.pt2.x * shrink_factor), this._mid_pt.y * (1 - shrink_factor) + (this.pt2.y * shrink_factor), null);
                    point_store pt3 = new point_store(-1, this._mid_pt.x * (1 - shrink_factor) + (this.pt3.x * shrink_factor), this._mid_pt.y * (1 - shrink_factor) + (this.pt3.y * shrink_factor), null);

                    this.shrunk_vertices[0] = pt1;
                    this.shrunk_vertices[1] = pt2;
                    this.shrunk_vertices[2] = pt3;
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

                public edge_store[] get_other_two_edges(edge_store the_edge)
                {
                    edge_store[] other_two_edge = new edge_store[2];

                    if (this.e1.Equals_without_orientation(the_edge) == true)
                    {
                        other_two_edge[0] = this.e2;
                        other_two_edge[1] = this.e3;
                        return other_two_edge;
                    }
                    else if (this.e2.Equals_without_orientation(the_edge) == true)
                    {
                        other_two_edge[0] = this.e3;
                        other_two_edge[1] = this.e1;
                        return other_two_edge;
                    }
                    else if (this.e3.Equals_without_orientation(the_edge) == true)
                    {
                        other_two_edge[0] = this.e1;
                        other_two_edge[1] = this.e2;
                        return other_two_edge;
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

                    if (has_same_sign(pab, pca) == false)
                    {
                        return false;
                    }


                    // Point test
                    if (e1.test_point_on_line(the_pt) == true)
                    {
                        return false;
                    }
                    else if (e2.test_point_on_line(the_pt) == true)
                    {
                        return false;
                    }
                    else if (e3.test_point_on_line(the_pt) == true)
                    {
                        return false;
                    }

                    return true;
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
