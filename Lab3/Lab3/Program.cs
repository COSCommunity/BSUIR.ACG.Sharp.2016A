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
                const int windowWidth = 800;
                const int windowHeight = 600;
                SDL.SDL_Init(SDL.SDL_INIT_EVERYTHING);
                IntPtr wnd = SDL.SDL_CreateWindow("Lab3", 100, 100, windowWidth, windowHeight,
                    SDL.SDL_WindowFlags.SDL_WINDOW_RESIZABLE |
                    SDL.SDL_WindowFlags.SDL_WINDOW_SHOWN);
                IntPtr renderer = SDL.SDL_CreateRenderer(wnd, -1, SDL.SDL_RendererFlags.SDL_RENDERER_ACCELERATED);

                var trapeze = new Polygon(
                    new[]
                    {
                        new SDL.SDL_Point() {x = 0, y = 0}, new SDL.SDL_Point() {x = 200, y = 0},
                        new SDL.SDL_Point() {x = 200, y = 200}, new SDL.SDL_Point() {x = 0, y = 100}
                    });
                var ellipse = new Ellipse(100, 200);
                var frame = new Polygon(
                    new[]
                    {
                        new SDL.SDL_Point() {x = 0, y = 0}, new SDL.SDL_Point() {x = windowWidth - 200, y = 0},
                        new SDL.SDL_Point() {x = windowWidth - 200, y = windowHeight - 200},
                        new SDL.SDL_Point() {x = 0, y = windowHeight - 200}
                    })
                {
                    TransformX = 100,
                    TransformY = 100,
                };
                DrawShapes(renderer, frame, trapeze, ellipse);

                SDL.SDL_Point trapezeVector = new SDL.SDL_Point() {x = -10, y = 10};
                SDL.SDL_Point ellipseVector = new SDL.SDL_Point() { x = -10, y = -5 };
                Random radnom = new Random();
                Timer timer = new Timer(e =>
                {
                    trapeze.TransformX += trapezeVector.x;
                    trapeze.TransformY += trapezeVector.y;
                    if (trapeze.TransformX >= windowWidth)
                    {
                        trapezeVector.x = -radnom.Next(0, 10);
                        trapezeVector.y = radnom.Next(0, 10);
                    }
                    else if (trapeze.TransformX <= 0)
                    {
                        trapezeVector.x = radnom.Next(0, 10);
                        trapezeVector.y = -radnom.Next(0, 10);
                    }
                    if (trapeze.TransformY >= windowHeight)
                    {
                        trapezeVector.x = radnom.Next(0, 10);
                        trapezeVector.y = -radnom.Next(0, 10);
                    }
                    else if (trapeze.TransformY <= 0)
                    {
                        trapezeVector.x = -radnom.Next(0, 10);
                        trapezeVector.y = radnom.Next(0, 10);
                    }
                    trapeze.Rotate += 10;

                    ellipse.TransformX += ellipseVector.x;
                    ellipse.TransformY += ellipseVector.y;
                    if (ellipse.TransformX >= windowWidth)
                    {
                        ellipseVector.x = -radnom.Next(0, 10);
                        ellipseVector.y = radnom.Next(0, 10);
                    }
                    else if (ellipse.TransformX <= 0)
                    {
                        ellipseVector.x = radnom.Next(0, 10);
                        ellipseVector.y = -radnom.Next(0, 10);
                    }
                    if (ellipse.TransformY >= windowHeight)
                    {
                        ellipseVector.x = -radnom.Next(0, 10);
                        ellipseVector.y = -radnom.Next(0, 10);
                    }
                    else if (ellipse.TransformY <= 0)
                    {
                        ellipseVector.x = radnom.Next(0, 10);
                        ellipseVector.y = radnom.Next(0, 10);
                    }
                    ellipse.Rotate += 5;

                    DrawShapes(renderer, frame, trapeze, ellipse);
                }, null, 0, 30);

                
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
                                    frame.TransformY += 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_UP:
                                    frame.TransformY -= 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_LEFT:
                                    frame.TransformX -= 10;
                                    break;
                                case SDL.SDL_Keycode.SDLK_RIGHT:
                                    frame.TransformX += 10;
                                    break;
                            }
                            DrawShapes(renderer, frame, trapeze, ellipse);
                            break;
                        }
                        case SDL.SDL_EventType.SDL_MOUSEBUTTONDOWN:
                        {
                            if (sdlEvent.button.button == SDL.SDL_BUTTON_LEFT)
                            {
                                if (frame.Scale < 1.2)
                                {
                                    frame.Scale += 0.05;
                                }
                            }
                            if (sdlEvent.button.button == SDL.SDL_BUTTON_RIGHT)
                            {
                                if (frame.Scale > 0.5)
                                {
                                    frame.Scale -= 0.05;
                                }
                            }
                            DrawShapes(renderer, frame, trapeze, ellipse);
                            break;
                        }

                    }                
                    Thread.Sleep(10);
                }
                SDL.SDL_DestroyRenderer(renderer);
                SDL.SDL_DestroyWindow(wnd);
                SDL.SDL_Quit();

            });
            thread.IsBackground = false;
            thread.Start();
        }

        private static object lockObject = new object();

        private static void DrawShapes(IntPtr renderer, Polygon frame, Polygon trapeze, Ellipse ellipse)
        {
            lock (lockObject)
            {
                SDL.SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
                SDL.SDL_RenderClear(renderer);
                SDL.SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
                DrawShape(renderer, frame);
                DrawShape(renderer, trapeze);
                DrawShape(renderer, ellipse);
                SDL.SDL_RenderPresent(renderer);
            }
        }

        private static void DrawShape(IntPtr renderer, Shape shape)
        {
            for (int i = 0; i < shape.ShapePoints.Length; i++)
            {
                DrawLine(renderer, shape.ShapePoints[i], shape.ShapePoints[(i + 1)%shape.ShapePoints.Length]);
            }
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
                set
                {
                    transformer.Rotate = value;
                    shapePoints = transformer.TransformPoints(Points);
                }
            }

            public double Scale
            {
                get { return transformer.Scale; }
                set
                {
                    transformer.Scale = value;
                    shapePoints = transformer.TransformPoints(Points);
                }
            }

            public int TransformX
            {
                get { return transformer.TransformX; }
                set
                {
                    transformer.TransformX = value;
                    shapePoints = transformer.TransformPoints(Points);
                }
            }

            public int TransformY
            {
                get { return transformer.TransformY; }
                set
                {
                    transformer.TransformY = value;
                    shapePoints = transformer.TransformPoints(Points);
                }
            }
            protected abstract SDL.SDL_Point[] Points { get; }

            private SDL.SDL_Point[] shapePoints;
            public SDL.SDL_Point[] ShapePoints => shapePoints ?? (shapePoints = transformer.TransformPoints(Points));
            public bool[] IsPointVisibleArray { get; protected set; }

            private readonly PointTransformer transformer = new PointTransformer()
            {
                Rotate = 0,
                Scale = 1,
                TransformX = 300,
                TransformY = 300
            };
        }

        private class Polygon : Shape
        {
            public Polygon(SDL.SDL_Point[] points)
            {
                this.Points = points;
                this.IsPointVisibleArray = points.Select(e => true).ToArray();
            }

            public Polygon(SDL.SDL_Point[] points, bool[] isPointVisibleArray)
            {
                this.Points = points;
                this.IsPointVisibleArray = isPointVisibleArray;
            }

            protected override SDL.SDL_Point[] Points { get; }
        }

        private class Ellipse : Shape, ICloneable
        {
            private const int PointNum = 10000;
            public int A { get; set; }
            public int B { get; set; }

            public Ellipse(int A, int B)
            {
                this.Points = InitPoints(A, B);
                this.IsPointVisibleArray = new bool[PointNum].Select(e => true).ToArray();
            }

            protected override SDL.SDL_Point[] Points { get; }

            private SDL.SDL_Point[] InitPoints(int A, int B)
            {

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

            public Ellipse Clone()
            {
                var el =  new Ellipse(A, B);
                el.Rotate = this.Rotate;
                el.Scale = this.Scale;
                el.TransformX = this.TransformX;
                el.TransformY = this.TransformY;
                return el;
            }

            object ICloneable.Clone()
            {
                return Clone();
            }
        }

        private static class IntersectionSearcher
        {
            public static Ellipse FindPolygonEllipseIntersection(Polygon polygon, Ellipse ellipse)
            {
                Ellipse result = ellipse.Clone();
                for (int i = 0; i < result.ShapePoints.Length; i++)
                {
                    if (IsPointInPolygon(polygon.ShapePoints, result.ShapePoints[i]))
                    {
                        result.IsPointVisibleArray[i] = false;
                    }
                }
                return result;
            }

            public static Polygon FindEllipsePolygonIntersection(Polygon polygon, Ellipse ellipse)
            {
                List<SDL.SDL_Point> interPoints = new List<SDL.SDL_Point>();
                List<bool> isPointVisible = new List<bool>();
                bool inPoly = true;
                for (int i = 0; i < ellipse.ShapePoints.Length; i++)
                {
                    if (inPoly)
                    {
                        if (IsPointInPolygon(polygon.ShapePoints, ellipse.ShapePoints[i]))
                        {

                            inPoly = false;
                        }
                        
                    }
                    else
                    {
                        if (!IsPointInPolygon(polygon.ShapePoints, ellipse.ShapePoints[i]))
                        {
                            inPoly = true;
                        }                        
                    }
                }
                return new Polygon(interPoints.ToArray(), isPointVisible.ToArray())
                {
                    Rotate = polygon.Rotate,
                    Scale = polygon.Scale,
                    TransformX = polygon.TransformX,
                    TransformY = polygon.TransformY
                };
            }

            private static bool IsPointInPolygon(SDL.SDL_Point[] polygonPoints, SDL.SDL_Point point)
            {
                bool result = false;
                for (int i = 0, j = polygonPoints.Length - 1; i < polygonPoints.Length; j = i++)
                {
                    if ((((polygonPoints[i].y <= point.y) && (point.y < polygonPoints[j].y)) ||
                         ((polygonPoints[j].y <= point.y) && (point.y < polygonPoints[i].y))) &&
                        (point.x >
                         (polygonPoints[j].x - polygonPoints[i].x)*(point.y - polygonPoints[i].y)/
                         (polygonPoints[j].y - polygonPoints[i].y) + polygonPoints[i].x))
                        result = !result;
                }
                return result;
            }

            public static Polygon FindPolygonsIntersection(Polygon lower, Polygon upper)
            {
                List<SDL.SDL_Point> interPoints = new List<SDL.SDL_Point>();
                List<bool> isPointVisible = new List<bool>();
                for (int i = 0; i < lower.ShapePoints.Length; i++)
                {
                    interPoints.Add(lower.ShapePoints[i]);
                    isPointVisible.Add(IsPointInPolygon(upper.ShapePoints, lower.ShapePoints[i]));
                    for (int j = 0; j < upper.ShapePoints.Length; j++)
                    {
                        var interPoint = FindLinesIntersectionPoint(lower.ShapePoints[i],
                            lower.ShapePoints[(i + 1)% lower.ShapePoints.Length], upper.ShapePoints[j],
                            upper.ShapePoints[(j + 1)% upper.ShapePoints.Length]);
                        if (interPoint.HasValue)
                        {
                            interPoints.Add(interPoint.Value);
                            isPointVisible.Add(false);
                        }
                    }
                }
                return new Polygon(interPoints.ToArray(), isPointVisible.ToArray())
                {
                    Rotate = lower.Rotate,
                    Scale = lower.Scale,
                    TransformX = lower.TransformX,
                    TransformY = lower.TransformY
                };
            }

            private static SDL.SDL_Point? FindLinesIntersectionPoint(SDL.SDL_Point f1, SDL.SDL_Point f2, SDL.SDL_Point s1, SDL.SDL_Point s2)
            {
                double k1 = 1.0*(f2.y - f1.y)/(f2.x - f2.y);
                double k2 = 1.0*(s2.y - s1.y)/(s2.x - s2.y);
                if (Math.Abs(k2 - k1) < 0.0000001)
                    return null;
                double c1 = 1.0*(f2.x*f1.y - f1.x*f2.y)/(f2.x - f2.y);
                double c2 = 1.0 * (s2.x * s1.y - s1.x * s2.y) / (s2.x - s2.y);
                var probX = (c1 - c2)/(k2 - k1);
                if (probX >= Math.Min(f1.x, f2.x) && probX <= Math.Max(f1.x, f2.x))
                {
                    return new SDL.SDL_Point() {x = (int) probX, y = (int) (k1*probX + c1)};
                }
                return null;
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
