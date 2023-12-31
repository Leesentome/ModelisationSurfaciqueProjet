#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>

void printLoadingBar(int progress, int total) {
    float percent = ((float)progress) / ((float)total);
    int filledBlocks = (int)(10 * percent);
    int emptyBlocks = 10 - filledBlocks;
    int digits = snprintf(NULL, 0, "%d", total);

    char progressBar[48 + 2 * digits];

    int length = snprintf(progressBar, sizeof(progressBar), "\r--> [%*d/%d] ", digits, progress, total);

    for (int i = 0; i < filledBlocks; i++) {
        length += snprintf(progressBar + length, sizeof(progressBar) - length, "█");
    }

    for (int i = 0; i < emptyBlocks; i++) {
        length += snprintf(progressBar + length, sizeof(progressBar) - length, "▒");
    }

    length += snprintf(progressBar + length, sizeof(progressBar) - length, " %6.2f%%", percent * 100);

    printf("%s", progressBar);
}

struct vertice {
    float x;
    float y;
    float z;
    int index;
    float Q[16];
};
int verticesLen;
int nbVertices;
struct vertice** vertices;

struct face {
    int len;
    int index;
    struct vertice** vertexIndex;
};
int facesLen;
int nbFaces;
struct face** faces;

int nbEdgesGiven;

void readFileOFF(const char* name) {
    FILE *file = fopen(name, "r");
    printf("File : \"%s\"\n", name);

    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char line[256];

    if (fgets(line, sizeof(line), file) != NULL) {
        printf(" └ Header: %s", line);
    } else {
        perror("Error reading header");
        exit(EXIT_FAILURE);
    }

    int numVertices, numFaces, numEdges;
    if (fscanf(file, "%d %d %d\n", &numVertices, &numFaces, &numEdges) != 3) {
        perror("Error reading counts");
        exit(EXIT_FAILURE);
    }

    printf("    └ Number of vertices: %d\n", numVertices);
    printf("    └ Number of faces: %d\n", numFaces);
    printf("    └ Number of edges: %d\n", numEdges);

    nbVertices = numVertices;
    verticesLen = numVertices;
    vertices = calloc(numVertices, sizeof(struct vertice*));

    for (int i = 0; i < numVertices; i++) {
        float x, y, z;
        if (fscanf(file, "%f %f %f\n", &x, &y, &z) != 3) {
            perror("Error reading vertex coordinates");
            exit(EXIT_FAILURE);
        }
        struct vertice* v = malloc(sizeof(struct vertice));
        v->x = x;
        v->y = y;
        v->z = z;
        v->index = i;
        for (int i = 0; i < 16; i++) v->Q[i] = 0;
        vertices[i] = v;
    }

    nbFaces = numFaces;
    facesLen = numFaces;
    faces = calloc(numFaces, sizeof(struct face*));

    for (int i = 0; i < numFaces; i++) {
        int numVerticesInFace;
        if (fscanf(file, "%d", &numVerticesInFace) != 1) {
            perror("Error reading number of vertices in face");
            exit(EXIT_FAILURE);
        }

        struct vertice** vertexIndexList = calloc(numVerticesInFace, sizeof(struct vertice*));

        for (int j = 0; j < numVerticesInFace; j++) {
            int vertexIndex;
            if (fscanf(file, "%d", &vertexIndex) != 1) {
                perror("Error reading vertex index in face");
                exit(EXIT_FAILURE);
            }
            vertexIndexList[j] = vertices[vertexIndex];
        }
        struct face* f = malloc(sizeof(struct face));
        f->len = numVerticesInFace;
        f->index = i;
        f->vertexIndex = vertexIndexList;
        faces[i] = f;
    }

    nbEdgesGiven = numEdges;

    fclose(file);
}

void writeFileOFF(const char* name) {
    FILE *file = fopen(name, "w");

    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "OFF\n");
    fprintf(file, "%i %i %i\n", nbVertices, nbFaces, 0);

    int index = 0;
    for (int i = 0; i < verticesLen; i++) {
        if (vertices[i] != NULL) {
            fprintf(file, "%f %f %f\n", vertices[i]->x, vertices[i]->y, vertices[i]->z);
            vertices[i]->index = index;
            index++;
        }
    }

    for (int i = 0; i < facesLen; i++) {
        if (faces[i] != NULL) {
            fprintf(file, "%i", faces[i]->len);

            for (int j = 0; j < faces[i]->len; j++) {
                fprintf(file, " %i", faces[i]->vertexIndex[j]->index);
            }

            fprintf(file, "\n");
        }
    }

    fclose(file);
}

struct halfEdge {
    struct vertice* startVertex;
    struct vertice* endVertex;
    struct face* face;
    struct halfEdge* prev;
    struct halfEdge* pair;
    struct halfEdge* next;
    float error;
    float vbar[3];
};

struct halfEdgeList {
    struct halfEdge* current;
    struct halfEdgeList* prev;
    struct halfEdgeList* next;
};
struct halfEdgeList* head = NULL;
struct halfEdgeList* edges = NULL;

void push(struct halfEdgeList** start, struct halfEdgeList* elt) {
    if (*start != NULL) {
        (*start)->prev = elt;
        elt->next = *start;
    }
    *start = elt;
}

void pop(struct halfEdgeList** start, struct halfEdgeList* elt) {
    if (elt == *start) {
        *start = elt->next;
    }
    if (elt->prev != NULL) {
        elt->prev->next = elt->next;
    }
    if (elt->next != NULL) {
        elt->next->prev = elt->prev;
    }
}

struct halfEdgeList* find(struct halfEdgeList* start, struct halfEdge* elt) {
    struct halfEdgeList* current = start;
    while (current != NULL) {
        if (current->current == elt) {
            return current;
        }
        current = current->next;
    }
    return NULL;
}

struct halfEdgeList* findPair(struct halfEdgeList* start, struct halfEdge* elt) {
    struct halfEdgeList* current = start;
    while (current != NULL) {

        if (elt->startVertex == current->current->endVertex &&
            elt->endVertex == current->current->startVertex) {
            return current;
        }

        current = current->next;
    }
    return NULL;
}

void computePlaneEquation(struct vertice* A, struct vertice* B, struct vertice* C, float* a, float* b, float* c, float* d) {
    float ABx = B->x - A->x;
    float ABy = B->y - A->y;
    float ABz = B->z - A->z;

    float ACx = C->x - A->x;
    float ACy = C->y - A->y;
    float ACz = C->z - A->z;

    *a = ABy * ACz - ABz * ACy;
    *b = ABz * ACx - ABx * ACz;
    *c = ABx * ACy - ABy * ACx;

    float length = sqrt(*a * *a + *b * *b + *c * *c);
    *a /= length;
    *b /= length;
    *c /= length;

    *d = -(*a * A->x + *b * A->y + *c * A->z);
}

void addQmatrix(struct vertice* vertex, float a, float b, float c, float d) {
    vertex->Q[0] += a * a;
    vertex->Q[1] += a * b;
    vertex->Q[2] += a * c;
    vertex->Q[3] += a * d;

    vertex->Q[4] += b * a;
    vertex->Q[5] += b * b;
    vertex->Q[6] += b * c;
    vertex->Q[7] += b * d;

    vertex->Q[8] += c * a;
    vertex->Q[9] += c * b;
    vertex->Q[10] += c * c;
    vertex->Q[11] += c * d;

    vertex->Q[12] += d * a;
    vertex->Q[13] += d * b;
    vertex->Q[14] += d * c;
    vertex->Q[15] += d * d;
}

void setHalfEdgeStruct() {
    printf("Creating and linking all halfEdge :\n");
    for (int i = 0; i < nbFaces; i++) {
        printLoadingBar(i, nbFaces);
        struct halfEdge* first = NULL;
        struct halfEdge* previous = NULL;
        struct halfEdge* current = NULL;

        struct vertice* A = faces[i]->vertexIndex[0];
        struct vertice* B = faces[i]->vertexIndex[1];
        struct vertice* C = faces[i]->vertexIndex[2];

        float a, b, c, d;
        computePlaneEquation(A, B, C, &a, &b, &c, &d);

        for (int j = 0; j < faces[i]->len; j++) {
            current = malloc(sizeof(struct halfEdge));
            current->startVertex = faces[i]->vertexIndex[j];
            current->endVertex = faces[i]->vertexIndex[(j+1)%faces[i]->len];
            current->face = faces[i];
            current->prev = previous;
            current->pair = NULL;
            current->next = NULL;
            current->error = 0;

            addQmatrix(faces[i]->vertexIndex[j], a, b, c, d);

            if (j == 0) {
                first = current;
            } else {
                previous->next = current;
            }

            struct halfEdgeList* nodePair = findPair(head, current);
            if (nodePair == NULL) {
                struct halfEdgeList* node = malloc(sizeof(struct halfEdgeList));
                node->current = current;
                node->prev=NULL;
                node->next=NULL;
                push(&head, node);
            } else {
                current->pair = nodePair->current;
                nodePair->current->pair = current;
                pop(&head, nodePair);
                free(nodePair);
            }

            // TOEDIT : Add both half edge -> maybe add only one

            struct halfEdgeList* edge = malloc(sizeof(struct halfEdgeList));
            edge->current = current;
            edge->prev=NULL;
            edge->next=NULL;
            push(&edges, edge);

            previous = current;
        }
        current->next = first;
        first->prev = current;
    }
    printLoadingBar(nbFaces, nbFaces);
    printf("\n\n");
}

float error(float* Q, float* v) {
    // Q*v1 = (a, b, c, d)
    float a = Q[ 0] * v[0] + Q[ 1] * v[1] + Q[ 2] * v[2]+ Q[ 3] * v[3];
    float b = Q[ 4] * v[0] + Q[ 5] * v[1] + Q[ 6] * v[2]+ Q[ 7] * v[3];
    float c = Q[ 8] * v[0] + Q[ 9] * v[1] + Q[10] * v[2]+ Q[11] * v[3];
    float d = Q[12] * v[0] + Q[13] * v[1] + Q[14] * v[2]+ Q[15] * v[3];

    return v[0] * a + v[1] * b + v[2] * c+ v[3] * d;
}

void computeHalfEdgeError(struct halfEdge* he) {
    float* Q1 = he->startVertex->Q;
    float* Q2 = he->endVertex->Q;

    float Qbar[16] = {0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0};
    
    for (int i = 0; i < 16; i++) Qbar[i] = Q1[i] + Q2[i];
        
    /*
        Qbar = a b c d     0  1  2  3
               e f g h     4  5  6  7
               i j k l     8  9 10 11
               0 0 0 1    12 13 14 15

        det = afk - agj - bek + bgi + cej - cfi
    */
    float det = (Qbar[0] * Qbar[5] * Qbar[10] + // + afk
                 Qbar[1] * Qbar[6] * Qbar[ 8] + // + bgi
                 Qbar[2] * Qbar[4] * Qbar[ 9] - // + cej
                 Qbar[2] * Qbar[5] * Qbar[ 8] - // - cfi
                 Qbar[1] * Qbar[4] * Qbar[10] - // - bek
                 Qbar[0] * Qbar[6] * Qbar[ 9]); // - agj
    
    if (det != 0) {
        /*
            Qbar-1 = α β ψ δ
                     ε φ γ η
                     ι ξ κ λ
                     0 0 0 1

            α = (fk - gj) / deno 
            β = (cj - bk) / deno
            ψ = (bg - cf) / deno
            δ = -(bgl - bhk - clf + chj + dfk - dgj) / deno
            ε = (gi - ek) / deno
            φ = (ak - ci) / deno
            γ = (ce - ag) / deno
            η = -(-agl + ahk + cel - chi - edk + dgi) / deno
            ι = (ej - fi) / deno
            ξ = (bi - aj) / deno
            κ = (af - be) / deno
            λ = -(afl - ahj - bel + bhi + dej - dfi) / deno

            μ = 1

            deno = afk - qgj - bek + bgi + cej - cfi
        */
        /*
            Qbar-1 * (0, 0, 0, 1) = (α β γ δ)

            α = -(bgl - bhk - clf + chj + dfk - dgj) / deno
            β =  (agl - ahk - cel + chi + edk - dgi) / deno
            γ = -(afl - ahj - bel + bhi + dej - dfi) / deno
            δ = 1

            deno = afk - qgj - bek + bgi + cej - cfi = det
        */
        float x_v_bar = -(Qbar[1] * Qbar[ 6] * Qbar[11] -        // + bgl
                          Qbar[1] * Qbar[ 7] * Qbar[10] -        // - bhk
                          Qbar[2] * Qbar[11] * Qbar[ 5] +        // - clf
                          Qbar[2] * Qbar[ 7] * Qbar[ 9] +        // + chj
                          Qbar[3] * Qbar[ 5] * Qbar[10] -        // + dfk
                          Qbar[3] * Qbar[ 6] * Qbar[ 9]) / det;  // - dgj
                            
        float y_v_bar =  (Qbar[0] * Qbar[ 6] * Qbar[11] -        // + agl
                          Qbar[0] * Qbar[ 7] * Qbar[10] -        // - ahk
                          Qbar[2] * Qbar[ 4] * Qbar[11] +        // - cel
                          Qbar[2] * Qbar[ 7] * Qbar[ 8] +        // + chi
                          Qbar[4] * Qbar[ 3] * Qbar[10] -        // + edk
                          Qbar[3] * Qbar[ 6] * Qbar[ 8]) / det;  // - dgi

        float z_v_bar = -(Qbar[0] * Qbar[ 5] * Qbar[11] -        // + afl
                          Qbar[0] * Qbar[ 7] * Qbar[ 9] -        // - ahj
                          Qbar[1] * Qbar[ 4] * Qbar[11] +        // - bel
                          Qbar[1] * Qbar[ 7] * Qbar[ 8] +        // + bhi
                          Qbar[3] * Qbar[ 4] * Qbar[ 9] -        // + dej
                          Qbar[3] * Qbar[ 5] * Qbar[ 8]) / det;  // - dfi
        he->vbar[0] = x_v_bar;
        he->vbar[1] = y_v_bar;
        he->vbar[2] = z_v_bar;

        float vbar[4] = {x_v_bar, y_v_bar, z_v_bar, 1};
        he->error = error(Qbar, vbar);
    } else {
        float v1[4] = {he->startVertex->x, he->startVertex->y, he->startVertex->z, 1};
        float v2[4] = {he->endVertex->x,   he->endVertex->y,   he->endVertex->z,   1};
        float v_mid[4] = {(v1[0]+v2[0])/2, (v1[1]+v2[1])/2, (v1[2]+v2[2])/2, 1};

        float err1 = error(Qbar, v1);
        float err2 = error(Qbar, v2);
        float err_mid = error(Qbar, v_mid);

        if (err1 < err2 && err1 < err_mid) {
            he->vbar[0] = v1[0];
            he->vbar[1] = v1[1];
            he->vbar[2] = v1[2];

            he->error = err1;
        } else if (err2 < err_mid) {
            he->vbar[0] = v2[0];
            he->vbar[1] = v2[1];
            he->vbar[2] = v2[2];

            he->error = err2;
        } else {
            he->vbar[0] = v_mid[0];
            he->vbar[1] = v_mid[1];
            he->vbar[2] = v_mid[2];

            he->error = err_mid;
        }
    }
}

void computeError() {
    struct halfEdgeList* current = edges;
    while (current != NULL) {

        computeHalfEdgeError(current->current);

        current = current->next;
    }
}

void clearGlobals() {
    for (int i = 0; i < verticesLen; i++) {
        if (vertices[i] != NULL) free(vertices[i]);
    }
    free(vertices);

    for (int i = 0; i < nbFaces; i++) {
        if (faces[i] != NULL) {
            free(faces[i]->vertexIndex);
            free(faces[i]);
        }
    }
    free(faces);

    while (edges != NULL) {
        struct halfEdgeList* current = edges;
        edges = edges->next;
        free(current);
    }
}

void contractEdge(struct halfEdge* he) {
    struct halfEdge* pair = he->pair;

    // suppression de la liste des aretes
    struct halfEdgeList* edge = find(edges, he);
    pop(&edges, edge);
    free(edge);
    struct halfEdgeList* pairEdge = find(edges, pair);
    pop(&edges, pairEdge);
    free(pairEdge);

    // ajustement des halfedge au point d'arrive
    struct halfEdge* cur = he->next;
    while (cur != pair) {
        cur->startVertex = he->startVertex;
        cur->pair->endVertex = pair->endVertex;
        cur = cur->pair->next;
        for (int j = 0; j < cur->face->len; j++) {
            if (cur->face->vertexIndex[j] == he->endVertex) {
                cur->face->vertexIndex[j] = he->startVertex;
            }  
        }
    }

    // ajustemet de la face directe
    struct face* heFace = he->face;
    if (heFace->len == 3) {
        he->prev->pair->pair = he->next->pair;
        he->next->pair->pair = he->prev->pair;

        struct halfEdgeList* prevEdge = find(edges, he->prev);
        pop(&edges, prevEdge);
        free(prevEdge);
        free(he->prev);
        struct halfEdgeList* nextEdge = find(edges, he->next);
        pop(&edges, nextEdge);
        free(nextEdge);
        free(he->next);

        nbFaces -= 1;
        faces[he->face->index] = NULL;
        free(he->face->vertexIndex);
        free(he->face);

    } else {
        he->prev->next = he->next;
        he->next->prev = he->prev;
        
        int j_ori = 0;
        int j_cop = 0;
        struct vertice** vertexIndexList = calloc(heFace->len - 1, sizeof(struct vertice*));
        while (j_ori < heFace->len) {
            if (heFace->vertexIndex[j_ori] != he->endVertex) {
                vertexIndexList[j_cop] = heFace->vertexIndex[j_ori];
                j_cop += 1;
            }
            j_ori += 1;
        }
        heFace->len--;
        free(heFace->vertexIndex);
        heFace->vertexIndex = vertexIndexList;
    }

    // ajustement de la face oposee
    struct face* pairFace = pair->face;
    if (pairFace->len == 3) {
        pair->prev->pair->pair = pair->next->pair;
        pair->next->pair->pair = pair->prev->pair;

        struct halfEdgeList* pairPrevEdge = find(edges, pair->prev);
        pop(&edges, pairPrevEdge);
        free(pairPrevEdge);
        free(pair->prev);
        struct halfEdgeList* pairNextEdge = find(edges, pair->next);
        pop(&edges, pairNextEdge);
        free(pairNextEdge);
        free(pair->next);
        
        nbFaces -= 1;
        faces[pair->face->index] = NULL;
        free(pair->face->vertexIndex);
        free(pair->face);

    } else {
        pair->prev->next = pair->next;
        pair->next->prev = pair->prev;

        int j_ori = 0;
        int j_cop = 0;
        struct vertice** vertexIndexList = calloc(pairFace->len - 1, sizeof(struct vertice*));
        while (j_ori < pairFace->len) {
            if (pairFace->vertexIndex[j_ori] != he->endVertex) {
                vertexIndexList[j_cop] = pairFace->vertexIndex[j_ori];
                j_cop += 1;
            }
            j_ori += 1;
        }
    }
    
    // suppression du sommet de la liste des sommets
    nbVertices -= 1;
    vertices[he->endVertex->index] = NULL;
    free(he->endVertex);

    free(pair);
    free(he);
}

float edgeLengthSquared(struct halfEdge* edge) {
    float x_len = edge->endVertex->x - edge->startVertex->x;
    float y_len = edge->endVertex->y - edge->startVertex->y;
    float z_len = edge->endVertex->z - edge->startVertex->z;
    return x_len * x_len + y_len * y_len + z_len * z_len;
}

struct halfEdge* findShortestEdge(struct halfEdgeList* start) {
    struct halfEdgeList* current = start;
    struct halfEdge* shortestEdge = NULL;
    while (current != NULL) {

        if (shortestEdge == NULL ||
            edgeLengthSquared(current->current) < edgeLengthSquared(shortestEdge)) {
            shortestEdge = current->current;
        }

        current = current->next;
    }
    return shortestEdge;
}

void contractShortestEdge() {
    struct halfEdge* shortestEdge = findShortestEdge(edges);
    contractEdge(shortestEdge);
}

struct halfEdge* findMinErrorEdge(struct halfEdgeList* start) {
    struct halfEdgeList* current = start;
    struct halfEdge* minErrorEdge = NULL;
    while (current != NULL) {

        if (minErrorEdge == NULL ||
            minErrorEdge->error > minErrorEdge->error) {
            minErrorEdge = current->current;
        }

        current = current->next;
    }
    return minErrorEdge;
}

void recomputeErrorVertex(struct halfEdge* sentinelEdgeForErrorRecompute) {
    struct halfEdge* current = sentinelEdgeForErrorRecompute->next->pair;
    computeHalfEdgeError(sentinelEdgeForErrorRecompute);
    sentinelEdgeForErrorRecompute->pair->error = sentinelEdgeForErrorRecompute->error;
    sentinelEdgeForErrorRecompute->pair->vbar[0] = sentinelEdgeForErrorRecompute->vbar[0];
    sentinelEdgeForErrorRecompute->pair->vbar[1] = sentinelEdgeForErrorRecompute->vbar[1];
    sentinelEdgeForErrorRecompute->pair->vbar[2] = sentinelEdgeForErrorRecompute->vbar[2];
    while (current != sentinelEdgeForErrorRecompute) {
        computeHalfEdgeError(current);
        current->pair->error = current->error;
        current->pair->vbar[0] = current->vbar[0];
        current->pair->vbar[1] = current->vbar[1];
        current->pair->vbar[2] = current->vbar[2];

        current = current->next->pair;
    }
}

void contractMinErrorEdge() {
    struct halfEdge* minErrorEdge = findMinErrorEdge(edges);
    struct halfEdge* sentinelEdgeForErrorRecompute = minErrorEdge->next->pair;
    contractEdge(minErrorEdge);
    recomputeErrorVertex(sentinelEdgeForErrorRecompute);
}

void contractEdgeTo(int goalVertices) {
    int nbIterations = nbVertices - goalVertices;
    printf("Contracting edge :\n");
    for (int i = 0; i < nbIterations; i++) {
        printLoadingBar(i, nbIterations);
        // contractShortestEdge();
        contractMinErrorEdge();
    }
    printLoadingBar(nbIterations, nbIterations);
    printf("\n\n");
}

int main(int argc, char *argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    // if (argc != 3) {
    //     // Print usage information if the number of arguments is incorrect
    //     printf("Usage: %s <input_file> <output_file>\n", argv[0]);
    //     return EXIT_FAILURE;
    // }

    // const char *fileFrom = argv[1];
    // const char *fileDest = argv[2];

    const char *fileFrom = "3D-Models/eight.off";
    const char *fileDest = "Simplified/eight_res_ter.off";

    readFileOFF(fileFrom);
    printf("\n");

    int goalVertices;
    printf("Enter the number of vertices for the output: ");
    scanf("%d", &goalVertices);
    printf("\n");

    setHalfEdgeStruct();

    computeError();
    contractEdgeTo(goalVertices);

    writeFileOFF(fileDest);

    clearGlobals();
    return 0;
}
