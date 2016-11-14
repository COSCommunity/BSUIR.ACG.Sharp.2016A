using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using SDL2;

namespace Lab2
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
                var shape = new Shape();
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
                                case SDL.SDL_Keycode.SDLK_9:
                                    shape.U = shape.U < 0.9 ? shape.U + 0.05 : shape.U;
                                    break;
                                case SDL.SDL_Keycode.SDLK_0:
                                    shape.U = shape.U > 0.1 ? shape.U - 0.05 : shape.U;
                                    break;
                                case SDL.SDL_Keycode.SDLK_1:
                                    shape.Angle = shape.Angle < 10 ? shape.Angle + 1 : shape.Angle;
                                    break;
                                case SDL.SDL_Keycode.SDLK_2:
                                    shape.Angle = shape.Angle > 3 ? shape.Angle - 1 : shape.Angle;
                                    break;
                                case SDL.SDL_Keycode.SDLK_TAB:
                                    shape.Mode = !shape.Mode;
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
        }

        private static void DrawShape(IntPtr renderer, Shape shape)
        {
            SDL.SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL.SDL_RenderClear(renderer);

            SDL.SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            var points = shape.Points;
            if (shape.Mode)
            {
                for (int j = 0; j < points.Length; j++)
                {
                    if (j%2 == 0)
                    {
                        var temp = points[j].ToList();
                        temp.Add(points[j][0]);
                        for (int i = 0; i < points[j].Length; i++)
                        {
                            DrawLine(renderer, temp[i], temp[i + 1]);
                        }
                    }
                    else
                    {
                        DrawCircle(renderer, points[j][0].x, points[j][0].y,
                            Convert.ToInt32(shape.Radiuses[j]*shape.Scale));
                    }
                    //SDL.SDL_RenderDrawLines(renderer, temp.ToArray(), arr.Length + 1);
                }
            }
            else
            {
                for (int j = 0; j < points.Length; j++)
                {
                    var temp = points[j].ToList();
                    temp.Add(points[j][0]);
                    for (int i = 0; i < points[j].Length; i++)
                    {
                        DrawLine(renderer, temp[i], temp[i + 1]);
                    }
                    //SDL.SDL_RenderDrawLines(renderer, temp.ToArray(), arr.Length + 1);
                }
            }

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

        private class Shape
        {
            private const int InnerShapeCount = 10;
            private int Radius = 100;

            private bool mode;

            public bool Mode
            {
                get { return mode; }
                set
                {
                    mode = value;
                    points = mode ? InitCircleShape() : InitShape();
                }
            }

            private int angle;

            public int Angle
            {
                get { return angle; }
                set
                {
                    angle = value;
                    points = InitShape();
                }
            }

            private double u;

            public double U
            {
                get { return u; }
                set
                {
                    u = value;
                    points = InitShape();
                }
            }

            public int Rotate { get; set; }
            public double Scale { get; set; }
            public int TransformX { get; set; }
            public int TransformY { get; set; }

            private SDL.SDL_Point[][] points;
            public SDL.SDL_Point[][] Points => TransformPoints();
            public int[] Radiuses = new int[InnerShapeCount];

            public Shape(int angles = 4)
            {
                u = 0.2;
                Scale = 1;
                Rotate = 0;
                TransformX = 300;
                TransformY = 300;
                angle = angles;
                points = InitShape();
            }

            private SDL.SDL_Point[][] InitCircleShape()
            {
                var result = new SDL.SDL_Point[InnerShapeCount][];
                var step = 2*Math.PI/4;
                double currentRadius = Radius;
                for (int i = 0; i < InnerShapeCount; i++)
                {
                    result[i] = new SDL.SDL_Point[4];
                    if (i%2 == 0)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            result[i][j].x = Convert.ToInt32(currentRadius*Math.Cos(step*j));
                            result[i][j].y = Convert.ToInt32(currentRadius*Math.Sin(step*j));
                        }
                        currentRadius = currentRadius/Math.Sqrt(2);
                    }
                    else
                    {
                        result[i][0].x = 0;
                        result[i][0].y = 0;
                        Radiuses[i] = Convert.ToInt32(currentRadius);
                    }
                }
                return result;
            }

            private SDL.SDL_Point[][] InitShape()
            {
                var result = new SDL.SDL_Point[InnerShapeCount][];
                var step = 2*Math.PI/Angle;
                for (int i = 0; i < InnerShapeCount; i++)
                {
                    result[i] = new SDL.SDL_Point[Angle];
                    for (int j = 0; j < Angle; j++)
                    {
                        if (i == 0)
                        {
                            result[i][j].x = Convert.ToInt32(Radius*Math.Cos(step*j));
                            result[i][j].y = Convert.ToInt32(Radius*Math.Sin(step*j));
                        }
                        else
                        {
                            result[i][j].x =
                                Convert.ToInt32((1 - U)*result[i - 1][j].x + U*result[i - 1][(j + 1)%Angle].x);
                            result[i][j].y =
                                Convert.ToInt32((1 - U)*result[i - 1][j].y + U*result[i - 1][(j + 1)%Angle].y);
                        }
                    }
                }
                return result;
            }

            private SDL.SDL_Point[][] TransformPoints()
            {
                var scale = Matrix.ScaleMatrix(Scale);
                var rotate = Matrix.RotateMatrix(Rotate);
                var transform = Matrix.TransformMatrix(TransformX, TransformY);
                var srtMatrix = transform*rotate*scale;
                var result = new SDL.SDL_Point[points.Length][];
                for (int i = 0; i < points.Length; i++)
                {
                    result[i] = new SDL.SDL_Point[points[i].Length];
                    for (int j = 0; j < result[i].Length; j++)
                    {
                        var vector = Matrix.VectorMatrix(points[i][j]);
                        result[i][j] = (SDL.SDL_Point) (srtMatrix*vector);
                    }
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
