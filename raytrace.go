package main

import (
	"os"
	"fmt"
	"math"
	"time"
	"math/rand"
	"image"
	"image/png"
	"image/color"
)

var BACKGROUND_COLOR = Color{135, 206, 235, 255}

type LightKind int
const (
	Ambient LightKind = iota
	Point
	Directional
)

type Canvas struct {
	width, height int
}

type Viewport struct {
	width, height, dist float64
}

type Vec3 struct {
	x, y, z float64
}

type Point3d = Vec3

type Direction = Vec3

type Color = color.RGBA

type Point2d struct {
	x, y float64
}

type Ray struct {
	origin Point3d
	dir Direction 
}

type Sphere struct {
	center Point3d
	radius float64
	color Color
	specular float64
	reflective float64
}

type Light struct {
	kind LightKind
	intensity float64
	position Point3d
	direction Direction
}

func CanvasToViewport(x, y int, c Canvas, v Viewport) Vec3 {
	vec := Vec3{float64(x) * v.width / float64(c.width), 
				float64(y) * v.height / float64(c.height),
				v.dist}
	return vec
}

func (a Vec3) Dot(b Vec3) float64 {
	return a.x * b.x + a.y * b.y + a.z * b.z
}

func (a Vec3) Add(b Vec3) Vec3 {
	return Vec3{a.x + b.x,
				a.y + b.y,
			    a.z + b.z}
}

func (a Vec3) Subtract(b Vec3) Vec3 {
	return Vec3{a.x - b.x,
				a.y - b.y,
			    a.z - b.z}
}

func (a Vec3) Scale(c float64) Vec3 {
	return Vec3{a.x * c,
				a.y * c,
			    a.z * c}
}

func (a Vec3) Multiply(b Vec3) Vec3 {
	return Vec3{a.x * b.x,
				a.y * b.y,
				a.z * b.z}
}

func (a Vec3) Magnitude() float64 {
	return math.Sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
}

func ScaleColor(color Color, c float64) Color {
	return Color{uint8(float64(color.R) * c),
			     uint8(float64(color.G) * c),
				 uint8(float64(color.B) * c),
				 color.A}
}

func AddColors(c1, c2 Color) Color {
	return Color{c1.R + c2.R,
				 c1.G + c2.G,
				 c1.B + c2.B,
				 c1.A}
}

func IntersectRaySphere(O Point3d, D Vec3, sphere Sphere) (float64, float64) {
	OC := O.Subtract(sphere.center)

	k1 := D.Dot(D)
	k2 := 2 * OC.Dot(D)
	k3 := OC.Dot(OC) - sphere.radius * sphere.radius

	discriminant := k2 * k2 - 4 * k1 * k3
	if discriminant < 0 {
		return math.Inf(1), math.Inf(1)
	}
	t1 := (-k2 + math.Sqrt(discriminant)) / (2 * k1)
	t2 := (-k2 - math.Sqrt(discriminant)) / (2 * k1)
	return t1, t2
}

func ClosestIntersection(O Point3d, D Vec3, tMin, tMax float64, scene *Scene) (Sphere, float64, bool) {
	minDist := math.Inf(1)
	var foundIntersect bool
	var closestSphere Sphere
	for _, sphere := range(scene.spheres) {
		t1, t2 := IntersectRaySphere(O, D, sphere)
		if tMin <= t1 && t1 < tMax {
			if t1 < minDist {
				minDist = t1
				closestSphere = sphere
				foundIntersect = true
			}
		}
		if tMin <= t2 && t2 < tMax {
			if t2 < minDist {
				minDist = t2
				closestSphere = sphere
				foundIntersect = true
			}
		}
	}
	return closestSphere, minDist, foundIntersect
}

func ComputeLighting(P Point3d, N, V Direction, s float64, scene *Scene) float64 {
	intensity := 0.0
	for _, light := range(scene.lights) {
		if light.kind == Ambient {
			intensity += light.intensity			
		} else {
			tMax := 1.0
			var L Direction
			if light.kind == Point {
				L = light.position.Subtract(P)
			}
			if light.kind == Directional {
				L = light.direction
				tMax = math.Inf(1)
			}

			// check for shadow
			eps := 1e-6 // machine epsilon is too small here
			_, _, isOccluded := ClosestIntersection(P, L, eps, tMax, scene)
			if isOccluded {
				continue
			}
		
			// Diffuse
			NL := N.Dot(L)
			if NL > 0 {
				intensity += light.intensity * NL / (N.Magnitude() * L.Magnitude())
			}
			//Specular
			if s != -1 {
				R := ReflectRay(L, N)
				RV := R.Dot(V)
				if RV > 0 {
					intensity += light.intensity * math.Pow(RV / (R.Magnitude() * V.Magnitude()), s)
				}
			}
		}
	}
	if intensity > 1 {
		return 1
	}
	return intensity
}

func TraceRay(O Point3d, D Vec3, tMin, tMax float64, scene *Scene, recursionDepth int) Color {
	closestSphere, dist, found := ClosestIntersection(O, D, tMin, tMax, scene)
	if !found {
		return BACKGROUND_COLOR
	}
	
	P := O.Add(D.Scale(dist))
	N := P.Subtract(closestSphere.center)
	N = N.Scale(1.0 / N.Magnitude())
	localColor := ScaleColor(closestSphere.color, ComputeLighting(P, N, D.Scale(-1), closestSphere.specular, scene))

	if recursionDepth <= 0 || closestSphere.reflective <= 0 {
		return localColor
	}

	R := ReflectRay(D.Scale(-1), N)
	eps := 1e-6 // machine espilon is too small here
	reflectedColor := TraceRay(P, R, eps, math.Inf(1), scene, recursionDepth - 1)

	return AddColors(ScaleColor(localColor, 1-closestSphere.reflective), 
	                 ScaleColor(reflectedColor, closestSphere.reflective))
}

func ReflectRay(R, N Vec3) Vec3 {
	NR := N.Dot(R)
	return N.Scale(2*NR).Subtract(R)
}

type Scene struct {
	spheres []Sphere
	lights []Light
}

func RandomColor() Color {
	r := uint8(rand.Intn(255))
	g := uint8(rand.Intn(255))
	b := uint8(rand.Intn(255))
	return Color{r, g, b, 255}
}

func RotateVectorX(theta float64, v Vec3) Vec3 {
	var R [3]Vec3
	R[0] = Vec3{1, 0, 0}
	R[1] = Vec3{0, math.Cos(theta), -math.Sin(theta)}
	R[2] = Vec3{0, math.Sin(theta), math.Cos(theta)}

	return Vec3{R[0].Dot(v),
				R[1].Dot(v),
				R[2].Dot(v)}
}

func RotateVectorY(theta float64, v Vec3) Vec3 {
	var R [3]Vec3
	R[0] = Vec3{math.Cos(theta), 0, math.Sin(theta)}
	R[1] = Vec3{0, 1, 0}
	R[2] = Vec3{-math.Sin(theta), 0, math.Cos(theta)}

	return Vec3{R[0].Dot(v),
				R[1].Dot(v),
				R[2].Dot(v)}
}

func RotateVectorZ(theta float64, v Vec3) Vec3 {
	var R [3]Vec3
	R[0] = Vec3{math.Cos(theta), -math.Sin(theta), 0}
	R[1] = Vec3{math.Sin(theta), math.Cos(theta), 0}
	R[2] = Vec3{0, 0, 1}

	return Vec3{R[0].Dot(v),
				R[1].Dot(v),
				R[2].Dot(v)}
}

func main() {
	fmt.Printf("Building Scene...\n")
	rand.Seed(time.Now().UnixNano())
	var spheres []Sphere
	for a := -11; a < 11; a++ {
		for b := 0; b < 11; b++ {
			center := Point3d{float64(a) + 0.9 * rand.Float64(), 0.2, float64(b) + 0.9 * rand.Float64()}
			color := RandomColor()
			choose_mat := rand.Float64()
			switch {
			case choose_mat > 0.95:
				// metal
				spheres = append(spheres, Sphere{center, 0.2, color, -1, 1})
			default:
				// diffuse
				specular := rand.Float64() * 1000
				reflective := rand.Float64() * 0.5
				spheres = append(spheres, Sphere{center, 0.2, color, specular, reflective})
			}
		}
	}

	spheres = append(spheres, Sphere{Point3d{0, -5000, 0}, 5000, Color{0, 154, 23, 255}, 1000, 0})

	l1 := Light{kind: Ambient, intensity: 0.3}
	l2 := Light{kind: Point, intensity: 0.7, position: Point3d{2, 1, 0}}
	l3 := Light{kind: Directional, intensity: 0.0, direction: Direction{1, 4, 1.3}}
	lights := []Light{l1, l2, l3}

	scene := Scene{spheres, lights}

	fmt.Printf("Ray Tracing!\n")
	lightX := -4.0
	lightPosList := []Point3d{}
	increment := 0.1
	for lightX <= 4 {
		lightPosList = append(lightPosList, Point3d{lightX, 2, 0})
		lightX = lightX + increment
	}

	
	for num, lightPos := range(lightPosList) {
		fmt.Printf("\r%d/%d", num, len(lightPosList))
		canv := Canvas{1920, 1080}

		upLeft := image.Point{0, 0}
		lowRight := image.Point{canv.width, canv.height}

		img := image.NewRGBA(image.Rectangle{upLeft, lowRight})

		origin := Point3d{0, 1, 0}
		
		view := Viewport{1, 1.0*float64(canv.height)/float64(canv.width), 0.75}
		reflectionDepth := 5
		scene.lights[1].position = lightPos
		for i := -canv.width / 2; i < canv.width / 2; i++ {
			for j := -canv.height / 2; j < canv.height / 2; j++ {
				D := CanvasToViewport(i, j, canv, view)
				D = RotateVectorX(math.Pi / 12, D)
				color := TraceRay(origin, D, 1, math.Inf(1), &scene, reflectionDepth)
				x := i + canv.width / 2
				y := -j + canv.height / 2
				img.Set(x, y, color)
			}
		}
		fname := fmt.Sprintf("image_%02d.png", num)
		f, _ := os.Create(fname)
		png.Encode(f, img)
	}
}