using System;
using System.Threading;
using System.Windows.Forms;
using SDL2;

namespace Lab1
{
    public partial class MainForm : Form
    {
        private IntPtr renderer;
        public MainForm()
        {
            InitializeComponent();
            Thread thread = new Thread(() =>
            {
                SDL.SDL_Init(SDL.SDL_INIT_EVERYTHING);
                IntPtr wnd = SDL.SDL_CreateWindow("Pascal shape SDL", 100, 100, 800, 600, SDL.SDL_WindowFlags.SDL_WINDOW_RESIZABLE |
                                                                                  SDL.SDL_WindowFlags.SDL_WINDOW_SHOWN);
                var shape = new PascalShape();
                renderer = SDL.SDL_CreateRenderer(wnd, -1, SDL.SDL_RendererFlags.SDL_RENDERER_ACCELERATED);
                DrawShape(shape);
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
                    DrawShape(shape);
                    Thread.Sleep(10);
                }
                SDL.SDL_DestroyRenderer(renderer);
                SDL.SDL_DestroyWindow(wnd);
                SDL.SDL_Quit();

            });
            thread.Start();
            thread.Join();
        }

        private void DrawShape(PascalShape shape)
        {
            SDL.SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL.SDL_RenderClear(renderer);
                   
            SDL.SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            var points = shape.Points;
            SDL.SDL_RenderDrawPoints(renderer, points, points.Length);
            SDL.SDL_RenderPresent(renderer);
        }   

        private class PascalShape
        {
            private const int PointArrayLength = 10000;
            private SDL.SDL_Point[] points;
            public int Rotate { get; set; }
            public double Scale { get; set; }
            public int TransformX { get; set; }
            public int TransformY { get; set; }            
            public SDL.SDL_Point[] Points => TransformPoints();

            public PascalShape()
            {
                Rotate = 0;
                Scale = 1;
                TransformX = 300;
                TransformY = 300;
                points = DrawPascal();
            }

            private SDL.SDL_Point[] TransformPoints()
            {
                var scale = Matrix.ScaleMatrix(Scale);
                var rotate = Matrix.RotateMatrix(Rotate);
                var transform = Matrix.TransformMatrix(TransformX, TransformY);
                var srtMatrix = transform*rotate*scale;
                var result = new SDL.SDL_Point[points.Length];
                for (int i = 0; i < result.Length; i++)
                {
                    var vector = Matrix.VectorMatrix(points[i]);
                    result[i] = (SDL.SDL_Point) (srtMatrix*vector);
                }
                return result;
            }

            private SDL.SDL_Point[] DrawPascal(int a = 100, int l = 50)
            {
                points = new SDL.SDL_Point[PointArrayLength];
                for (int i = 0; i < PointArrayLength; i++)
                {
                    var t = 2*Math.PI/PointArrayLength*i;
                    double x = a*Math.Pow(Math.Cos(t), 2) + l*Math.Cos(t);
                    double y = a*Math.Cos(t)*Math.Sin(t) + l*Math.Sin(t);
                    points[i] = new SDL.SDL_Point()
                    {
                        x = Convert.ToInt32(x),
                        y = Convert.ToInt32(y)
                    };
                }
                return points;
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
                double piDegree = degree / 180.0 * Math.PI;
                return new Matrix()
                {
                    Values =
                        new[]
                        {
                            new[] {Math.Cos(piDegree), Math.Sin(piDegree), 0}, new[] {-Math.Sin(piDegree), Math.Cos(piDegree), 0},
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
                    Values = new[] { new double[] { point.x, 0, 0 }, new double[] { point.y, 0, 0 }, new double[] { 1, 0, 0 } }
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

        private void MainForm_Shown(object sender, EventArgs e)
        {
            Hide();
            Close();
        }
    }
}
