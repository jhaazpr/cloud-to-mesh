#include "CGL/CGL.h"

#include "collada.h"
#include "meshEdit.h"
#include "bezierPatch.h"
#include "mergeVertices.h"
#include "shaderUtils.h"

#include <iostream>

using namespace std;
using namespace CGL;

#define msg(s) cerr << "[Collada Viewer] " << s << endl;

bool loadPLY(const char * path, std::vector<Vector3D> & out_vertices,
              std::vector<Polygon> &out_polygons);

int loadFile(MeshEdit* collada_viewer, const char* path) {

    Scene* scene = new Scene();

    std::string path_str = path;
    if (path_str.substr(path_str.length() - 4, 4) == ".dae") {
        if (ColladaParser::load(path, scene) < 0) {
            delete scene;
            return -1;
        }
    } else if (path_str.substr(path_str.length() - 4, 4) == ".bez") {
        Camera* cam = new Camera();
        cam->type = CAMERA;
        Node node;
        node.instance = cam;
        scene->nodes.push_back(node);
        Polymesh* mesh = new Polymesh();

        FILE* file = fopen(path, "r");
        int n = 0;
        fscanf(file, "%d", &n);
        for (int i = 0; i < n; i++) {
            BezierPatch patch;
            patch.loadControlPoints(file);
            patch.add2mesh(mesh);
            mergeVertices(mesh);
        }
        fclose(file);

        mesh->type = POLYMESH;
        node.instance = mesh;
        scene->nodes.push_back(node);
    } else if (path_str.substr(path_str.length() - 4, 4) == ".ply") {
        //NOTE: Jasper's code begins here
        cout << "You cool" << endl;
        // Use GL_POINTS primitve somehow
        std::vector<Vector3D> out_vertices;
        std::vector<Polygon> out_polygons;
        bool success = loadPLY(path, out_vertices, out_polygons);
        //TODO: actually use the points

        // Set up scene
        Camera* cam = new Camera();
        cam->type = CAMERA;
        Node node;
        node.instance = cam;
        scene->nodes.push_back(node);
        PointCloud* point_cloud = new PointCloud();
        Polymesh *polymesh = new Polymesh();

        // Add vertices to scene
        for (Vector3D v : out_vertices) {
          point_cloud->vertices.push_back(v);
          polymesh->vertices.push_back(v);
        }

        // Add polygons to scene
        for (Polygon p : out_polygons) {
          polymesh->polygons.push_back(p);
        }

        point_cloud->type = POINT_CLOUD;
        node.instance = point_cloud;
        scene->nodes.push_back(node);

        // polymesh->type = POLYMESH;
        // node.instance = polymesh;
        // scene->nodes.push_back(node);

    } else {
        return -1;
    }

    collada_viewer->load(scene);

    GLuint tex = makeTex("envmap/envmap.png");
    if (!tex) tex = makeTex("../envmap/envmap.png");
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, tex);
    glActiveTexture(GL_TEXTURE2);

    return 0;
}

bool loadPLY(const char * path, std::vector<Vector3D> & out_vertices,
                std::vector<Polygon> &out_polygons) {
    printf("Load PLY called\n");
    FILE* file = fopen(path, "r");
    std::vector<Vector3D> temp_vertices;
    std::vector<Polygon> temp_polygons;
    int vertex_count = 0;
    int face_count = 0;

    // Read until we get to element vertices <number>
    char lineHeader[512];
    while (true) {

        // read the first word of the line
        int res = fscanf(file, "%s", lineHeader);
//        printf("res: %s\n", lineHeader);
        if (res == EOF) {
//            printf("Hit EOF\n");
            break; // EOF = End Of File. Quit the loop.
        }

        if (strcmp(lineHeader, "end_header") == 0) {
            break;
        }

        if ( strcmp( lineHeader, "element" ) == 0 ) {
            char type[128];
            int count;
            fscanf(file, "%s %d\n", type, &count );
            if (strcmp(type, "vertex") == 0) {
                vertex_count = count;
            }

            if (strcmp(type, "face") == 0) {
                face_count = count;
            }
        }

    }

    cout << "LoadPLY: read number of vertices: " << vertex_count << endl;
    // cout << "LoadPLY: read number of faces: " << face_count << endl;

    // Main loop of parsing vertices
    while (--vertex_count >= 0) {
        Vector3D vertex;
        float dummy0, dummy1;

        //NOTE: uncomment one based on the format of the ply file
        // fscanf(file, "%lf %lf %lf %f %f\n", &vertex.x, &vertex.y, &vertex.z , &dummy0, &dummy1);
        fscanf(file, "%lf %lf %lf\n", &vertex.x, &vertex.y, &vertex.z);

        temp_vertices.push_back(vertex);
    }

    // // Main loop of parsing faces
    // while (--face_count >= 0) {
    //     Polygon polygon;
    //     Index idx0, idx1, idx2;
    //     float dummy0;
    //
    //     //NOTE: uncomment one based on the format of the ply file
    //     // fscanf(file, "%lf %lf %lf %f %f\n", &vertex.x, &vertex.y, &vertex.z , &dummy0, &dummy1);
    //     fscanf(file, "%f %lu %lu %lu \n", &dummy0, &idx0, &idx1, &idx2);
    //     polygon.vertex_indices.push_back(idx0);
    //     polygon.vertex_indices.push_back(idx1);
    //     polygon.vertex_indices.push_back(idx2);
    //     temp_polygons.push_back(polygon);
    // }

    out_vertices = temp_vertices;
    out_polygons = temp_polygons;
    return true;
}

int main(int argc, char** argv) {

    // create viewer
    Viewer viewer = Viewer();

    // create collada_viewer
    MeshEdit* collada_viewer = new MeshEdit();

    // set collada_viewer as renderer
    viewer.set_renderer(collada_viewer);

    // init viewer
    viewer.init();

    // load tests
    if (argc == 2) {
        if (loadFile(collada_viewer, argv[1]) < 0) exit(0);
    } else {
        msg("Usage: collada_viewer <path to scene file>");
        exit(0);
    }

    // start viewer
    viewer.start();

    return 0;
}
