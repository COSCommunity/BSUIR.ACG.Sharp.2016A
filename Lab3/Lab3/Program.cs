using System;
using System.Collections;
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
                    }) {Layer = 1};
                var ellipse = new Ellipse(100, 200) {Layer = 2};
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
                    Layer = -5
                };

                List<Shape> shapes = new List<Shape>() {frame, trapeze, ellipse};

                DrawShapes(renderer, shapes);

                SDL.SDL_Point trapezeVector = new SDL.SDL_Point() {x = -10, y = 10};
                SDL.SDL_Point ellipseVector = new SDL.SDL_Point() { x = -10, y = -10 };
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
                    trapeze.Rotate += 5;

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
                    ellipse.Rotate += 3;
          
                    DrawShapes(renderer, shapes);
                }, null, 0, 30);
                
                bool quit = false;
                bool stopped = false;
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
                                case SDL.SDL_Keycode.SDLK_SPACE:
                                    if (!stopped)
                                    {
                                        timer.Change(0, Timeout.Infinite);
                                        stopped = true;
                                    }
                                    else
                                    {
                                        timer.Change(0, 30);
                                        stopped = false;
                                    }
                                    break;
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
                            break;
                        }

                    }                
                    Thread.Sleep(30);
                }
                SDL.SDL_DestroyRenderer(renderer);
                SDL.SDL_DestroyWindow(wnd);
                SDL.SDL_Quit();

            });
            thread.IsBackground = false;
            thread.Start();
        }

        private static object lockObject = new object();

        private static void DrawShapes(IntPtr renderer, IEnumerable<Shape> shapes)
        {
            lock (lockObject)
            {
                SDL.SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
                SDL.SDL_RenderClear(renderer);
                SDL.SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
                var interShapes = IntersectionSearcher.FindShapesIntersction(shapes);
                foreach (var shape in interShapes)
                    DrawShape(renderer, shape);
                SDL.SDL_RenderPresent(renderer);
            }
        }

        private static void DrawShape(IntPtr renderer, InterceptedShape shape)
        {
            for (int i = 0; i < shape.Points.Length; i++)
            {
                if (!shape.Points[i].IsVisible && !shape.Points[(i + 1) % shape.Points.Length].IsVisible)
                {
                    if (!(shape.Points.Length > 50 && i % 2 == 0))
                        DrawDashLine(renderer, shape.Points[i].Point, shape.Points[(i + 1) % shape.Points.Length].Point);
                }
                else
                    DrawLine(renderer, shape.Points[i].Point, shape.Points[(i + 1) % shape.Points.Length].Point);
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

        private static void DrawDashLine(IntPtr renderer, SDL.SDL_Point p1, SDL.SDL_Point p2)
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
                    int count = 0;
                    for (int x = p1.x + sx, y = p1.y, i = 1; i <= dx; i++, x += sx)
                    {
                        if (d > 0)
                        {
                            d += d2;
                            y += sy;
                        }
                        else
                            d += d1;
                        if ((count/10)%2 == 0)
                            SDL.SDL_RenderDrawPoint(renderer, x, y);
                        count++;
                    }
                }
                else
                {
                    int d = (dx << 1) - dy;
                    int d1 = dx << 1;
                    int d2 = (dx - dy) << 1;
                    SDL.SDL_RenderDrawPoint(renderer, p1.x, p1.y);
                    int count = 0;
                    for (int x = p1.x, y = p1.y + sy, i = 1; i <= dy; i++, y += sy)
                    {
                        if (d > 0)
                        {
                            d += d2;
                            x += sx;
                        }
                        else
                            d += d1;
                        if ((count/10)%2 == 0)
                            SDL.SDL_RenderDrawPoint(renderer, x, y);
                        count++;
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

        private class VisibleSdlPoint
        {
            public SDL.SDL_Point Point { get; set; }
            public bool IsVisible { get; set; }
        }


        private static class IntersectionSearcher
        {

            public static IEnumerable<InterceptedShape> FindShapesIntersction(IEnumerable<Shape> shapes)
            {
                var temp = shapes as Shape[] ?? shapes.ToArray();
                return temp.Select(e => FindIntersectionedShape(e, temp));
            }

            public static InterceptedShape FindIntersectionedShape(Shape shape, IEnumerable<Shape> shapes)
            {
                var interShape =
                    new InterceptedShape(
                        shape.Points.Select(e => new VisibleSdlPoint() { Point = e, IsVisible = true }).ToArray());
                foreach (var temp in shapes.Where(e => Math.Abs(e.Layer) > Math.Abs(shape.Layer)))
                {
                    interShape = FindShapeIntersection(interShape, temp);
                }
                return interShape;
            }

            public static InterceptedShape FindShapeIntersection(InterceptedShape lower, Shape upper)
            {
                List<VisibleSdlPoint> interPoints = new List<VisibleSdlPoint>();
                for (int i = 0; i < lower.Points.Length; i++)
                {
                    bool r = !IsPointInPolygon(upper.Points, lower.Points[i].Point);
                    if (upper.Layer < 0)
                        r = !r;
                    interPoints.Add(new VisibleSdlPoint()
                    {
                        Point = lower.Points[i].Point,
                        IsVisible = lower.Points[i].IsVisible && r
                    });
                    for (int j = 0; j < upper.Points.Length; j++)
                    {
                        var interPoint = FindLinesIntersectionPoint(lower.Points[i].Point,
                            lower.Points[(i + 1) % lower.Points.Length].Point, upper.Points[j],
                            upper.Points[(j + 1) % upper.Points.Length]);
                        if (interPoint.HasValue)
                        {
                            interPoints.Add(new VisibleSdlPoint()
                            {
                                Point = interPoint.Value,
                                IsVisible = false
                            });
                        }
                    }
                }
                return new InterceptedShape(interPoints.ToArray());
            }

            private static SDL.SDL_Point? FindLinesIntersectionPoint(SDL.SDL_Point f1, SDL.SDL_Point f2, SDL.SDL_Point s1, SDL.SDL_Point s2)
            {
                double k1 = 1.0*(f2.y - f1.y)/(f2.x - f1.x);
                double k2 = 1.0*(s2.y - s1.y)/(s2.x - s1.x);
                if (Math.Abs(k2 - k1) < 0.0000001)
                    return null;
                double c1 = f1.y - k1 * f1.x;
                double c2 = s1.y - k2 * s1.x;
                var probX = (c2 - c1)/(k1 - k2);
                if (probX >= Math.Min(f1.x,f2.x) && probX <= Math.Max(f1.x,f2.x) && probX >= Math.Min(s1.x,s2.x) && probX <= Math.Max(s1.x,f2.x))
                {
                    return new SDL.SDL_Point() {x = Convert.ToInt32(probX), y = Convert.ToInt32((k1*probX + c1))};
                }
                return null;
            }

            private static bool IsPointInPolygon(SDL.SDL_Point[] polygonPoints, SDL.SDL_Point point)
            {
                bool result = false;
                for (int i = 0, j = polygonPoints.Length - 1; i < polygonPoints.Length; j = i++)
                {
                    if ((((polygonPoints[i].y <= point.y) && (point.y < polygonPoints[j].y)) ||
                         ((polygonPoints[j].y <= point.y) && (point.y < polygonPoints[i].y))) &&
                        (point.x >
                         (polygonPoints[j].x - polygonPoints[i].x) * (point.y - polygonPoints[i].y) /
                         (polygonPoints[j].y - polygonPoints[i].y) + polygonPoints[i].x))
                        result = !result;
                }
                return result;
            }

            #region unusefulMethods
            //public static Ellipse FindPolygonEllipseIntersection(Polygon polygon, Ellipse ellipse)
            //{
            //    Ellipse result = ellipse.Clone();
            //    for (int i = 0; i < result.ShapePoints.Length; i++)
            //    {
            //        if (IsPointInPolygon(polygon.ShapePoints, result.ShapePoints[i]))
            //        {
            //            result.IsPointVisibleArray[i] = false;
            //        }
            //    }
            //    return result;
            //}

            //public static Polygon FindEllipsePolygonIntersection(Polygon polygon, Ellipse ellipse)
            //{
            //    List<SDL.SDL_Point> interPoints = new List<SDL.SDL_Point>();               
            //    bool inPoly = true;
            //    for (int i = 0; i < ellipse.ShapePoints.Length; i++)
            //    {
            //        if (inPoly)
            //        {
            //            if (IsPointInPolygon(polygon.ShapePoints, ellipse.ShapePoints[i]))
            //            {
            //                interPoints.Add(ellipse.ShapePoints[i]);
            //                inPoly = false;
            //            }

            //        }
            //        else
            //        {
            //            if (!IsPointInPolygon(polygon.ShapePoints, ellipse.ShapePoints[i]))
            //            {
            //                interPoints.Add(ellipse.ShapePoints[i - 1]);
            //                inPoly = true;
            //            }
            //        }
            //    }
            //    List<SDL.SDL_Point> joinedPoints = new List<SDL.SDL_Point>();
            //    List<bool> isPointVisible = new List<bool>();
            //    for (int i = 0; i < polygon.ShapePoints.Length; i++)
            //    {
            //        joinedPoints.Add(polygon.ShapePoints[i]);
            //        isPointVisible.Add(!IsPointInPolygon(ellipse.ShapePoints, polygon.ShapePoints[i]));

            //        foreach (var p in interPoints)
            //        {
            //            var index = 0;
            //            for (int j = 1; j < polygon.ShapePoints.Length; j++)
            //            {
            //                if (Length(p, polygon.ShapePoints[j]) < Length(p, polygon.ShapePoints[index]))
            //                {
            //                    index = j;
            //                }
            //            }
            //            if (index == i)
            //            {
            //                joinedPoints.Add(p);
            //                isPointVisible.Add(false);
            //            }
            //        }
            //    }
            //    return new Polygon(joinedPoints.ToArray(), isPointVisible.ToArray())
            //    {
            //        Rotate = polygon.Rotate,
            //        Scale = polygon.Scale,
            //        TransformX = polygon.TransformX,
            //        TransformY = polygon.TransformY
            //    };
            //}

            //private static double Length(SDL.SDL_Point p1, SDL.SDL_Point p2)
            //{
            //    return Math.Sqrt(Math.Pow(p1.x - p2.x, 2) + Math.Pow(p1.y - p2.y, 2));
            //}
            #endregion
        }

        private abstract class Shape
        {
            private SDL.SDL_Point[] shapePoints;
            protected SDL.SDL_Point[] initPoints;

            private readonly PointTransformer transformer = new PointTransformer()
            {
                Rotate = 0,
                Scale = 1,
                TransformX = 0,
                TransformY = 0
            };

            public int Rotate
            {
                get { return transformer.Rotate; }
                set
                {
                    transformer.Rotate = value;
                    shapePoints = transformer.TransformPoints(initPoints);
                }
            }

            public double Scale
            {
                get { return transformer.Scale; }
                set
                {
                    transformer.Scale = value;
                    shapePoints = transformer.TransformPoints(initPoints);
                }
            }

            public int TransformX
            {
                get { return transformer.TransformX; }
                set
                {
                    transformer.TransformX = value;
                    shapePoints = transformer.TransformPoints(initPoints);
                }
            }

            public int TransformY
            {
                get { return transformer.TransformY; }
                set
                {
                    transformer.TransformY = value;
                    shapePoints = transformer.TransformPoints(initPoints);
                }
            }

            public int Layer { get; set; }

            public virtual SDL.SDL_Point[] Points
                => shapePoints ?? (shapePoints = transformer.TransformPoints(initPoints));
        }

        private class InterceptedShape
        {
            public InterceptedShape(VisibleSdlPoint[] points)
            {
                this.Points = points;
            }
            public VisibleSdlPoint[] Points { get; }
            public int Layer { get; set; }
        }

        private class Polygon : Shape
        {
            public Polygon(SDL.SDL_Point[] points)
            {
                this.initPoints = points;
            }
        }

        private class Ellipse : Shape
        {
            private const int PointNum = 100;

            public Ellipse(int A, int B)
            {
                this.initPoints = InitPoints(A, B);
            }

            private SDL.SDL_Point[] InitPoints(int A, int B)
            {

                var result = new SDL.SDL_Point[PointNum];
                for (int i = 0; i < result.Length; i++)
                {
                    result[i] = new SDL.SDL_Point
                    {
                        x = Convert.ToInt32(A * Math.Cos(2 * Math.PI / PointNum * i)),
                        y = Convert.ToInt32(B * Math.Sin(2 * Math.PI / PointNum * i))
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
