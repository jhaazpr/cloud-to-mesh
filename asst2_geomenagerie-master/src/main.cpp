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

bool loadPLY(const char * path, std::vector<Vector3D> & out_vertices);

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
        bool success = loadPLY(path, out_vertices);
        //TODO: actually use the points

        // Set up scene
        Camera* cam = new Camera();
        cam->type = CAMERA;
        Node node;
        node.instance = cam;
        scene->nodes.push_back(node);
        PointCloud* point_cloud = new PointCloud();

        // Add vertices to scene
        for (Vector3D v : out_vertices) {
          point_cloud->vertices.push_back(v);
        }

        point_cloud->type = POINT_CLOUD;
        node.instance = point_cloud;
        scene->nodes.push_back(node);
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

bool loadPLY(const char * path, std::vector<Vector3D> & out_vertices) {
    printf("Load PLY called\n");
    FILE* file = fopen(path, "r");
    std::vector<Vector3D> temp_vertices;
    int vertex_count = 0;

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

        if ( strcmp( lineHeader, "element" ) == 0 ) {
            char type[128];
            int count;
            fscanf(file, "%s %d\n", type, &count );
            if (strcmp(type, "vertex") == 0) {
                vertex_count = count;
            }
        }

        if (strcmp(lineHeader, "end_header") == 0) {
            break;
        }
    }

    cout << "LoadPLY: read number of vertices: " << vertex_count << endl;

    // Main loop of parsing vertices
    while (vertex_count-- > 0) {
        Vector3D vertex;
        float dummy0, dummy1;
        fscanf(file, "%lf %lf %lf %f %f\n", &vertex.x, &vertex.y, &vertex.z , &dummy0, &dummy1);
        temp_vertices.push_back(vertex);
//        printf("%f %f %f\n", vertex.x, vertex.y, vertex.z);
    }

    out_vertices = temp_vertices;
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
