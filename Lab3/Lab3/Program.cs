using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using SDL2;

namespace Lab3
{
    class Program
    {
        static void Main(string[] args)
        {
            Thread thread = new Thread(() =>
            {
                SDL.SDL_Init(SDL.SDL_INIT_EVERYTHING);
                IntPtr wnd = SDL.SDL_CreateWindow("Pascal shape SDL", 100, 100, 800, 600,
                    SDL.SDL_WindowFlags.SDL_WINDOW_RESIZABLE |
                    SDL.SDL_WindowFlags.SDL_WINDOW_SHOWN);

                IntPtr renderer = SDL.SDL_CreateRenderer(wnd, -1, SDL.SDL_RendererFlags.SDL_RENDERER_ACCELERATED);
                var shape = new Polygon(1, new SDL.SDL_Point[] {});
                DrawShape(renderer, shape);
                bool quit = false;
                while (!quit)
                {
                    SDL.SDL_Event sdlEvent;
                    SDL.SDL_PollEvent(out sdlEvent);
                    switch (sdlEvent.type)
                    {
                        case SDL.SDL_EventType.SDL_QUIT:
                        {
                            quit = true;
                            break;
                        }
                        case SDL.SDL_EventType.SDL_KEYDOWN:
                        {
                            var key = sdlEvent.key;
                            switch (key.keysym.sym)
                            {
                                case SDL.SDL_Keycode.SDLK_DOWN:
                                    shape.TransformY += 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_UP:
                                    shape.TransformY -= 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_LEFT:
                                    shape.TransformX -= 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_RIGHT:
                                    shape.TransformX += 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_LCTRL:
                                    shape.Rotate += 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_RCTRL:
                                    shape.Rotate -= 10;
                                    break;
                            }
                            break;
                        }
                        case SDL.SDL_EventType.SDL_MOUSEBUTTONDOWN:
                        {
                            if (sdlEvent.button.button == SDL.SDL_BUTTON_LEFT)
                            {
                                if (shape.Scale < 3)
                                {
                                    shape.Scale += 0.1;
                                }
                            }
                            if (sdlEvent.button.button == SDL.SDL_BUTTON_RIGHT)
                            {
                                if (shape.Scale > 0.2)
                                {
                                    shape.Scale -= 0.1;
                                }
                            }
                            break;
                        }

                    }
                    DrawShape(renderer, shape);
                    Thread.Sleep(10);
                }
                SDL.SDL_DestroyRenderer(renderer);
                SDL.SDL_DestroyWindow(wnd);
                SDL.SDL_Quit();
            });
            thread.IsBackground = false;
            thread.Start();
            //
        }

        private static void DrawShape(IntPtr renderer, Shape shape)
        {
            SDL.SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL.SDL_RenderClear(renderer);
            SDL.SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);



            SDL.SDL_RenderPresent(renderer);
        }

        private static void DrawLine(IntPtr renderer, SDL.SDL_Point p1, SDL.SDL_Point p2)
        {
            {
                int dx = Math.Abs(p2.x - p1.x);
                int dy = Math.Abs(p2.y - p1.y);
                int sx = p2.x >= p1.x ? 1 : -1;
                int sy = p2.y >= p1.y ? 1 : -1;
                if (dy <= dx)
                {
                    int d = (dy << 1) - dx;
                    int d1 = dy << 1;
                    int d2 = (dy - dx) << 1;
                    SDL.SDL_RenderDrawPoint(renderer, p1.x, p1.y);
                    for (int x = p1.x + sx, y = p1.y, i = 1; i <= dx; i++, x += sx)
                    {
                        if (d > 0)
                        {
                            d += d2;
                            y += sy;
                        }
                        else
                            d += d1;
                        SDL.SDL_RenderDrawPoint(renderer, x, y);
                    }
                }
                else
                {
                    int d = (dx << 1) - dy;
                    int d1 = dx << 1;
                    int d2 = (dx - dy) << 1;
                    SDL.SDL_RenderDrawPoint(renderer, p1.x, p1.y);
                    for (int x = p1.x, y = p1.y + sy, i = 1; i <= dy; i++, y += sy)
                    {
                        if (d > 0)
                        {
                            d += d2;
                            x += sx;
                        }
                        else
                            d += d1;
                        SDL.SDL_RenderDrawPoint(renderer, x, y);
                    }
                }
            }
        }

        private static void DrawCircle(IntPtr renderer, int _x, int _y, int radius)
        {
            int x = 0, y = radius, gap = 0, delta = (2 - 2*radius);
            while (y >= 0)
            {
                SDL.SDL_RenderDrawPoint(renderer, _x + x, _y + y);
                SDL.SDL_RenderDrawPoint(renderer, _x + x, _y - y);
                SDL.SDL_RenderDrawPoint(renderer, _x - x, _y - y);
                SDL.SDL_RenderDrawPoint(renderer, _x - x, _y + y);
                gap = 2*(delta + y) - 1;
                if (delta < 0 && gap <= 0)
                {
                    x++;
                    delta += 2*x + 1;
                    continue;
                }
                if (delta > 0 && gap > 0)
                {
                    y--;
                    delta -= 2*y + 1;
                    continue;
                }
                x++;
                delta += 2*(x - y);
                y--;
            }
        }

        private abstract class Shape
        {
            public int Rotate
            {
                get { return transformer.Rotate; }
                set { transformer.Rotate = value; }
            }

            public double Scale
            {
                get { return transformer.Scale; }
                set { transformer.Scale = value; }
            }

            public int TransformX
            {
                get { return transformer.TransformX; }
                set { transformer.TransformX = value; }
            }

            public int TransformY
            {
                get { return transformer.TransformY; }
                set { transformer.TransformY = value; }
            }

            public int Layer { get; set; }
            protected abstract SDL.SDL_Point[] Points { get; }

            private readonly PointTransformer transformer = new PointTransformer()
            {
                Rotate = 0,
                Scale = 1,
                TransformX = 300,
                TransformY = 300
            };

            public SDL.SDL_Point[] GetShapePoints()
            {
                return transformer.TransformPoints(Points);
            }
        }

        private class Polygon : Shape
        {
            public Polygon(int layer, SDL.SDL_Point[] points)
            {
                this.Layer = layer;
                this.Points = points;
            }

            protected override SDL.SDL_Point[] Points { get; }
        }

        private class Ellipse : Shape
        {
            public Ellipse(int layer, int A, int B)
            {
                this.Layer = layer;
                this.Points = InitPoints(A, B);
            }

            protected override SDL.SDL_Point[] Points { get; }

            private SDL.SDL_Point[] InitPoints(int A, int B)
            {
                const int PointNum = 10000;
                var result = new SDL.SDL_Point[PointNum];
                for (int i = 0; i < result.Length; i++)
                {
                    result[i] = new SDL.SDL_Point
                    {
                        x = Convert.ToInt32(A*Math.Cos(2*Math.PI/PointNum*i)),
                        y = Convert.ToInt32(B*Math.Sin(2*Math.PI/PointNum*i))
                    };
                }
                return result;
            }
        }

        private class PointTransformer
        {
            public double Scale { get; set; }
            public int Rotate { get; set; }
            public int TransformX { get; set; }
            public int TransformY { get; set; }

            public SDL.SDL_Point[] TransformPoints(SDL.SDL_Point[] points)
            {
                var scale = Matrix.ScaleMatrix(Scale);
                var rotate = Matrix.RotateMatrix(Rotate);
                var transform = Matrix.TransformMatrix(TransformX, TransformY);
                var srtMatrix = transform*rotate*scale;
                var result = new SDL.SDL_Point[points.Length];
                for (int i = 0; i < points.Length; i++)
                {
                    var vector = Matrix.VectorMatrix(points[i]);
                    result[i] = (SDL.SDL_Point) (srtMatrix*vector);
                }
                return result;
            }
        }

        private class Matrix
        {
            private double[][] Values { get; set; }

            public static Matrix ScaleMatrix(double pow)
            {

                return new Matrix()
                {
                    Values = new[] {new[] {pow, 0, 0}, new[] {0, pow, 0}, new double[] {0, 0, 1}}
                };
            }

            public static Matrix RotateMatrix(int degree)
            {
                double piDegree = degree/180.0*Math.PI;
                return new Matrix()
                {
                    Values =
                        new[]
                        {
                            new[] {Math.Cos(piDegree), Math.Sin(piDegree), 0},
                            new[] {-Math.Sin(piDegree), Math.Cos(piDegree), 0},
                            new double[] {0, 0, 1}
                        }
                };
            }

            public static Matrix TransformMatrix(int x, int y)
            {
                return new Matrix()
                {
                    Values = new[] {new double[] {1, 0, x}, new double[] {0, 1, y}, new double[] {0, 0, 1}}
                };
            }

            public static Matrix VectorMatrix(SDL.SDL_Point point)
            {
                return new Matrix()
                {
                    Values = new[] {new double[] {point.x, 0, 0}, new double[] {point.y, 0, 0}, new double[] {1, 0, 0}}
                };
            }

            public static Matrix operator *(Matrix matrix1, Matrix matrix2)
            {
                Matrix result = new Matrix() {Values = new[] {new double[3], new double[3], new double[3]}};
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        result.Values[i][j] = 0;
                        for (int k = 0; k < 3; k++)
                        {
                            result.Values[i][j] += matrix1.Values[i][k]*matrix2.Values[k][j];
                        }
                    }
                }
                return result;
            }

            public static explicit operator SDL.SDL_Point(Matrix matrix)
            {
                return new SDL.SDL_Point()
                {
                    x = Convert.ToInt32(matrix.Values[0][0]),
                    y = Convert.ToInt32(matrix.Values[1][0])
                };
            }
        }
    }
}
