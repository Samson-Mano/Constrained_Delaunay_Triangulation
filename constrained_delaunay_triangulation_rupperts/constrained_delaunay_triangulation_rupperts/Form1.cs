using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace constrained_delaunay_triangulation
{
    public partial class Form1 : Form
    {
        map_generator main_map = new map_generator();
        pslg_datastructure pslg_data = new pslg_datastructure();

        public Form1()
        {
            InitializeComponent();
        }

        private void button_randommap_Click(object sender, EventArgs e)
        {
            the_static_class.bounding_box = new SizeF((float)(0.8 * main_pic.Width), (float)(0.8 * main_pic.Height));
            the_static_class.bounding_midpt = new PointF((float)(0.5 * main_pic.Width), (float)(0.5 * main_pic.Height));
            generate_map();
        }

        public void generate_map()
        {
            // generate map 
            main_map = new map_generator(the_static_class.bounding_box.Width, the_static_class.bounding_box.Height);
            main_map.detect_surfaces(ref pslg_data.all_surfaces);
            mt_pic.Refresh();
        }

        public static class the_static_class
        {
            public static SizeF bounding_box;
            public static PointF bounding_midpt;

            public static bool ispaint_label = false; // static variable to control whether to paint id or not

            public static bool is_animate_checked = false; // static variable to control animation timing;
            public static bool is_paint_incircle = false;
            public static bool is_paint_mesh = true;

            public static int instance_counter_at; // static variable to control instances count
            public static int inpt_timer_interval = 500; // 0.5 seconds, static variable to control the interval of the timer

            public static double B_var = Math.Sqrt(2); // Ruppert's algorithm condition B_var = Math.Sqrt(2) && Paul chew's second algorithm condition (Math.Sqrt(5) * 0.5)

            /// <summary>
                        /// Function to Check the valid of Numerical text from textbox.text
                        /// </summary>
                        /// <param name="tB_txt">Textbox.text value</param>
                        /// <param name="Negative_check">Is negative number Not allowed (True) or allowed (False)</param>
                        /// <param name="zero_check">Is zero Not allowed (True) or allowed (False)</param>
                        /// <returns>Return the validity (True means its valid) </returns>
                        /// <remarks></remarks>
            public static bool test_a_textboxvalue_validity_int(string tb_txt, bool n_chk, bool z_chk)
            {
                bool is_valid = false;
                //This function returns false if the textbox doesn't contains number 
                if (Int32.TryParse(tb_txt, out Int32 number) == true)
                {
                    is_valid = true;

                    if (n_chk == true) // check for negative number
                    {
                        if (Convert.ToInt32(tb_txt) < 0)
                        {
                            is_valid = false;
                        }
                    }

                    if (z_chk == true) // check for zero number
                    {
                        if (Convert.ToInt32(tb_txt) == 0)
                        {
                            is_valid = false;
                        }
                    }
                }
                return is_valid;
            }

            /// <summary>
                        /// Function to convert double to single (mostly used in System.Drawing functions)
                        /// </summary>
                        /// <param name="value"></param>
                        /// <returns></returns>
            public static float to_single(double value)
            {
                return (float)value;
            }

            /// <summary>
                        /// Funtion to check NAN or Infinity value
                        /// </summary>
                        /// <param name="chkval"></param>
                        /// <returns></returns>
            public static bool Isval_NAN_or_Infinity(double chkval)
            {
                return (double.IsNaN(chkval) || double.IsInfinity(chkval));
            }

            public static Color GetRandomColor(int hash)
            {
                //Random randonGen = new Random(DateTime.Now.Millisecond.GetHashCode());
                Random randomGen = new Random((hash+19)* DateTime.Now.Millisecond.GetHashCode());
                Color randomColor = Color.FromArgb(randomGen.Next(0,256), randomGen.Next(0,256), randomGen.Next(0,256));
                return randomColor;
            }

            public static System.Drawing.Drawing2D.HatchStyle GetRandomHatchStyle(int hash)
            {
                //Random randomGen = new Random(hash.GetHashCode());
                //int randomhatchindex = randomGen.Next(0, 6);
                int randomhatchindex = hash > 6 ? 0 : hash;
                System.Drawing.Drawing2D.HatchStyle style = new System.Drawing.Drawing2D.HatchStyle();

                switch (randomhatchindex)
                {
                    case 0:
                        style = System.Drawing.Drawing2D.HatchStyle.BackwardDiagonal;
                        break;
                    case 1:
                        style = System.Drawing.Drawing2D.HatchStyle.DashedVertical;
                        break;
                    case 2:
                        style = System.Drawing.Drawing2D.HatchStyle.Cross;
                        break;
                    case 3:
                        style = System.Drawing.Drawing2D.HatchStyle.DiagonalCross;
                        break;
                    case 4:
                        style = System.Drawing.Drawing2D.HatchStyle.HorizontalBrick;
                        break;
                    case 5:
                        style = System.Drawing.Drawing2D.HatchStyle.LightDownwardDiagonal;
                        break;
                    case 6:
                        style = System.Drawing.Drawing2D.HatchStyle.LightUpwardDiagonal;
                        break;
                    default:
                        break;
                }

                return style;
            }

            // Return the angle ABC.
            // Return a value between PI and -PI.
            // Note that the value is the opposite of what you might
            // expect because Y coordinates increase downward.
            public static double GetAngle(double Ax, double Ay,
                double Bx, double By, double Cx, double Cy)
            {
                // Get the dot product.
                double dot_product = DotProduct(Ax, Ay, Bx, By, Cx, Cy);

                // Get the cross product.
                double cross_product = CrossProductLength(Ax, Ay, Bx, By, Cx, Cy);

                // Calculate the angle.
                return Math.Atan2(cross_product, dot_product);
            }

            // Return the cross product AB x BC.
            // The cross product is a vector perpendicular to AB
            // and BC having length |AB| * |BC| * Sin(theta) and
            // with direction given by the right-hand rule.
            // For two vectors in the X-Y plane, the result is a
            // vector with X and Y components 0 so the Z component
            // gives the vector's length and direction.
            public static double CrossProductLength(double Ax, double Ay,
                double Bx, double By, double Cx, double Cy)
            {
                // Get the vectors' coordinates.
                double BAx = Ax - Bx;
                double BAy = Ay - By;
                double BCx = Cx - Bx;
                double BCy = Cy - By;

                // Calculate the Z coordinate of the cross product.
                return (BAx * BCy - BAy * BCx);
            }

            // Return the dot product AB · BC.
            // Note that AB · BC = |AB| * |BC| * Cos(theta).
            private static double DotProduct(double Ax, double Ay,
                double Bx, double By, double Cx, double Cy)
            {
                // Get the vectors' coordinates.
                double BAx = Ax - Bx;
                double BAy = Ay - By;
                double BCx = Cx - Bx;
                double BCy = Cy - By;

                // Calculate the dot product.
                return (BAx * BCx + BAy * BCy);
            }

            /// <summary>
            /// Determines if the given point is inside the polygon
            /// </summary>
            /// <param name="polygon">the vertices of polygon</param>
            /// <param name="testPoint">the given point</param>
            /// <returns>true if the point is inside the polygon; otherwise, false</returns>
            public static bool IsPointInPolygon_1(PointF[] polygon, PointF testPoint)
            {
                bool result = false;
                int j = polygon.Count() - 1;
                for (int i = 0; i < polygon.Count(); i++)
                {
                    if (polygon[i].Y < testPoint.Y && polygon[j].Y >= testPoint.Y || polygon[j].Y < testPoint.Y && polygon[i].Y >= testPoint.Y)
                    {
                        if (polygon[i].X + (testPoint.Y - polygon[i].Y) / (polygon[j].Y - polygon[i].Y) * (polygon[j].X - polygon[i].X) < testPoint.X)
                        {
                            result = !result;
                        }
                    }
                    j = i;
                }
                return result;
            }

            // Return True if the point is in the polygon.
            public static bool IsPointInPolygon(PointF[] polygon, PointF testPoint)
            {
                // Get the angle between the point and the
                // first and last vertices.
                int max_point = polygon.Length - 1;
                double total_angle = GetAngle(
                    polygon[max_point].X, polygon[max_point].Y,
                    testPoint.X, testPoint.Y,
                    polygon[0].X, polygon[0].Y);

                // Add the angles from the point
                // to each other pair of vertices.
                for (int i = 0; i < max_point; i++)
                {
                    total_angle += GetAngle(
                        polygon[i].X, polygon[i].Y,
                        testPoint.X, testPoint.Y,
                        polygon[i + 1].X, polygon[i + 1].Y);
                }

                // The total angle should be 2 * PI or -2 * PI if
                // the point is in the polygon and close to zero
                // if the point is outside the polygon.
                // The following statement was changed. See the comments.
                //return (Math.Abs(total_angle) > 0.000001);
                return (Math.Abs(total_angle) > 1);
            }


            public static double SignedPolygonArea(PointF[] polygon)
            {
                //The total calculated area is negative if the polygon is oriented clockwise

                // Add the first point to the end.
                int num_points = polygon.Length;
                PointF[] pts = new PointF[num_points + 1];
                polygon.CopyTo(pts, 0);
                pts[num_points] = polygon[0];

                // Get the areas.
                double area = 0;
                for (int i = 0; i < num_points; i++)
                {
                    area +=
                        (pts[i + 1].X - pts[i].X) *
                        (pts[i + 1].Y + pts[i].Y) / 2;
                }

                // Return the result.
                return area;
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            the_static_class.bounding_box = new SizeF((float)(0.5 * main_pic.Width), (float)(0.5 * main_pic.Height));
            the_static_class.bounding_midpt = new PointF((float)(0.5 * main_pic.Width), (float)(0.5 * main_pic.Height));
            //generate_map();
        }

        private void Form1_SizeChanged(object sender, EventArgs e)
        {
            the_static_class.bounding_box = new SizeF((float)(0.8 * main_pic.Width), (float)(0.8 * main_pic.Height));
            the_static_class.bounding_midpt = new PointF((float)(0.5 * main_pic.Width), (float)(0.5 * main_pic.Height));
            mt_pic.Refresh();
        }

        private void main_pic_Paint(object sender, PaintEventArgs e)
        {
            e.Graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality; // Paint high quality
            e.Graphics.TranslateTransform(the_static_class.bounding_midpt.X, the_static_class.bounding_midpt.Y); // Translate transform to make the orgin as center

            //e.Graphics.DrawRectangle(new Pen(Color.Brown, 2), 0, 0, 8, 8);
            System.Drawing.Graphics gr1 = e.Graphics;
            //main_map.paint_me(ref gr1);
            pslg_data.paint_me(ref gr1);
        }

        private void main_pic_MouseClick(object sender, MouseEventArgs e)
        {
            // mouse click event to create or delete mesh
            if (e.Button == MouseButtons.Left)
            {
                // create the mesh
                pslg_data.select_and_set_mesh(e.X - the_static_class.bounding_midpt.X, e.Y - the_static_class.bounding_midpt.Y, true);
            }

            if (e.Button == MouseButtons.Right)
            {
                // delete the mesh
                pslg_data.select_and_set_mesh(e.X - the_static_class.bounding_midpt.X, e.Y - the_static_class.bounding_midpt.Y, false);
            }

            mt_pic.Refresh();
        }

        private void main_pic_MouseMove(object sender, MouseEventArgs e)
        {
            label_xylocation.Text = "(" + (e.X- the_static_class.bounding_midpt.X).ToString() + ", " + (e.Y- the_static_class.bounding_midpt.Y).ToString() + " )";
            label_xylocation.Refresh();
        }
    }
}
