#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

int maxi(int a, int b) {
    return (a > b) ? a : b;
}

void printLoadingBar(int progress, int total) {
    float percent = ((float)progress) / ((float)total);
    int filledBlocks = (int)(10 * percent);
    int emptyBlocks = 10 - filledBlocks;
    int digits = snprintf(NULL, 0, "%d", total);

    // char progressBar[48 + 2 * digits];
    int totLength = 48 + 2 * digits;
    char *progressBar = (char *)malloc(totLength*sizeof(char));

    int length = snprintf(progressBar, totLength, "\r--> [%*d/%d] ", digits, progress, total);

    for (int i = 0; i < filledBlocks; i++) {
        length += snprintf(progressBar + length, totLength - length, "█");
    }

    for (int i = 0; i < emptyBlocks; i++) {
        length += snprintf(progressBar + length, totLength - length, "▒");
    }

    length += snprintf(progressBar + length, totLength - length, " %6.2f%%", percent * 100);

    progressBar[length] = '\0';

    printf("%s", progressBar);
    free(progressBar);
}

struct vertice {
    float x;
    float y;
    float z;
    int index;
    float* Q;
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
    printf("Reading file : \"%s\"\n", name);

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
        v->Q = calloc(16, sizeof(float));
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
    printf("Writing file : \"%s\"\n", name);
    printf(" └ Number of vertices: %d\n", nbVertices);
    printf(" └ Number of faces: %d\n", nbFaces);

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
    printf("Nb vertex : %d\n", index);

    index = 0;
    for (int i = 0; i < facesLen; i++) {
        if (faces[i] != NULL) {
            fprintf(file, "%i", faces[i]->len);

            for (int j = 0; j < faces[i]->len; j++) {
                fprintf(file, " %i", faces[i]->vertexIndex[j]->index);
            }

            fprintf(file, "\n");
            index++;
        }
    }
    printf("Nb face : %d\n", index);


    fclose(file);
}

struct halfEdgeError {
    bool valid;
    bool good;
    float err;
};

bool compareErrors(struct halfEdgeError* a, struct halfEdgeError* b) {
    // return a->err < b->err;
    if (a->valid && a->good) {
        if (b->valid && b->good) {
            return a->err < b->err;
        } else {
            return true;
        }
    } else if (a->valid) {
        if (b->valid && b->good) {
            return false;
        } else if (b->valid) {
            return a->err < b->err;
        } else {
            return true;
        }
    } else if (a->good) {
        if (b->valid) {
            return false;
        } else if (b->good) {
            return a->err < b->err;
        } else {
            return true;
        }
    } else {
        if (b->valid || b->good) {
            return false;
        } else {
            return a->err < b->err;
        }
    }
}

struct halfEdge {
    struct vertice* startVertex;
    struct vertice* endVertex;
    struct face* face;
    struct halfEdge* prev;
    struct halfEdge* pair;
    struct halfEdge* next;
    struct halfEdgeError* error;
    float* vbar;
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

struct halfEdgeAVL {
    struct halfEdge* data;
    struct halfEdgeAVL* parent;
    struct halfEdgeAVL* equal;
    struct halfEdgeAVL* left;
    struct halfEdgeAVL* right;
    int height;
};
struct halfEdgeAVL* edgesTree;

int height(struct halfEdgeAVL* node) {
    return (node == NULL) ? 0 : node->height;
}

struct halfEdgeAVL* newNode(struct halfEdge* hEdge) {
    struct halfEdgeAVL* node = malloc(sizeof(struct halfEdgeAVL));
    node->data = hEdge;
    node->parent = NULL;
    node->equal = NULL;
    node->left = NULL;
    node->right = NULL;
    node->height = 1;
    return node;
}

struct halfEdgeAVL* rightRotate(struct halfEdgeAVL* tree) {
    struct halfEdgeAVL* lSubTree = tree->left;
    struct halfEdgeAVL* T2 = lSubTree->right;

    lSubTree->right = tree;
    lSubTree->parent = tree->parent;
    tree->left = T2;
    tree->parent = lSubTree;

    tree->height = maxi(height(tree->left), height(tree->right)) + 1;
    lSubTree->height = maxi(height(lSubTree->left), height(lSubTree->right)) + 1;

    return lSubTree;
}

struct halfEdgeAVL* leftRotate(struct halfEdgeAVL* tree) {
    struct halfEdgeAVL* rSubTree = tree->right;
    struct halfEdgeAVL* T2 = rSubTree->left;

    rSubTree->left = tree;
    rSubTree->parent = tree->parent;
    tree->right = T2;
    tree->parent = rSubTree;

    tree->height = maxi(height(tree->left), height(tree->right)) + 1;
    rSubTree->height = maxi(height(rSubTree->left), height(rSubTree->right)) + 1;

    return rSubTree;
}

int getBalance(struct halfEdgeAVL* tree) {
    return (tree == NULL) ? 0 : height(tree->left) - height(tree->right);
}

struct halfEdgeAVL* minValueNode(struct halfEdgeAVL* tree) {
    struct halfEdgeAVL* current = tree;

    while (current->left != NULL) {
        current = current->left;
    }

    return current;
}

struct halfEdgeAVL* insertNode(struct halfEdgeAVL* tree, struct halfEdge* hEdge) {
    if (tree == NULL) {
        return newNode(hEdge);
    }

    if (compareErrors(hEdge->error, tree->data->error)) {
        tree->left = insertNode(tree->left, hEdge);
        tree->left->parent = tree;
    } else if (compareErrors(tree->data->error, hEdge->error)) {
        tree->right = insertNode(tree->right, hEdge);
        tree->right->parent = tree;
    } else {
        struct halfEdgeAVL* current = tree;
        while (current->equal != NULL) {
            current = current->equal;
        }
        current->equal = newNode(hEdge);

        return tree;
    }

    tree->height = 1 + maxi(height(tree->left), height(tree->right));

    int balance = getBalance(tree);

    // Left Left Case
    if (balance > 1 && compareErrors(hEdge->error, tree->left->data->error)) {
        return rightRotate(tree);
    }

    // Right Right Case
    if (balance < -1 && compareErrors(tree->right->data->error, hEdge->error)) {
        return leftRotate(tree);
    }

    // Left Right Case
    if (balance > 1 && compareErrors(tree->left->data->error, hEdge->error)) {
        tree->left = leftRotate(tree->left);
        return rightRotate(tree);
    }

    // Right Left Case
    if (balance < -1 && compareErrors(hEdge->error, tree->right->data->error)) {
        tree->right = rightRotate(tree->right);
        return leftRotate(tree);
    }

    return tree;
}

struct halfEdgeAVL* findNode(struct halfEdgeAVL* tree, struct halfEdge* hEdge) {
    while (tree != NULL) {
        if (tree->data == hEdge) {
            return tree;
        } else if (compareErrors(hEdge->error, tree->data->error)) {
            tree = tree->left;
        } else if (compareErrors(tree->data->error, hEdge->error)) {
            tree = tree->right;
        } else {
            struct halfEdgeAVL* current = tree->equal;
            while (current != NULL) {
                if (current->data == hEdge) {
                    return current;
                }
                current = current->equal;
            }
            return NULL;
        }
    }
    return NULL;
}

struct halfEdgeAVL* deleteNode(struct halfEdgeAVL* tree, struct halfEdgeAVL* node) {
    if (tree == node) {
        if (tree->equal) {
            tree->data = tree->equal->data;
            struct halfEdgeAVL* temp = tree->equal;
            tree->equal = tree->equal->equal;
            free(temp);
            return tree;
        } else if (tree->left == NULL || tree->right == NULL) {
            struct halfEdgeAVL* temp = tree->left ? tree->left : tree->right;

            if (temp == NULL) {
                temp = tree;
                tree = NULL;
            } else {
                tree->data = temp->data;
                tree->parent = temp->parent;
                tree->equal = temp->equal;
                tree->left = temp->left;
                tree->right = temp->right;
                tree->height = temp->height;
            }

            free(temp);
        } else {
            struct halfEdgeAVL* temp = minValueNode(tree->right);
            struct halfEdge* treeData = tree->data;
            tree->data = temp->data;
            tree->equal = temp->equal;
            
            temp->equal = NULL;
            temp->data = treeData;

            tree->right = deleteNode(tree->right, temp);
        }
    } else if (compareErrors(node->data->error, tree->data->error)) {
        tree->left = deleteNode(tree->left, node);
    } else if (compareErrors(tree->data->error, node->data->error)) {
        tree->right = deleteNode(tree->right, node);
    } else {
        struct halfEdgeAVL* current = tree;
        while (current->equal != NULL) {
            if (current->equal == node) {
                current->equal = current->equal->equal;
                free(node);
                break;
            }
            current = current->equal;
        }
        return tree;
    }

    if (tree == NULL) {
        return tree;
    }

    tree->height = 1 + maxi(height(tree->left), height(tree->right));

    int balance = getBalance(tree);

    // Left Left Case
    if (balance > 1 && getBalance(tree->left) >= 0) {
        return rightRotate(tree);
    }

    // Left Right Case
    if (balance > 1 && getBalance(tree->left) < 0) {
        tree->left = leftRotate(tree->left);
        return rightRotate(tree);
    }

    // Right Right Case
    if (balance < -1 && getBalance(tree->right) <= 0) {
        return leftRotate(tree);
    }

    // Right Left Case
    if (balance < -1 && getBalance(tree->right) > 0) {
        tree->right = rightRotate(tree->right);
        return leftRotate(tree);
    }

    return tree;
}

void updateTreeFromNode(struct halfEdgeAVL* modified, struct halfEdgeError* oldErr) {
    struct halfEdge* nodeData = modified->data;
    struct halfEdgeError* newErr = modified->data->error;
    nodeData->error = oldErr;
 
    edgesTree = deleteNode(edgesTree, modified);

    nodeData->error = newErr;

    edgesTree = insertNode(edgesTree, nodeData);

}

void freeAVLTree(struct halfEdgeAVL* tree) {
    if (tree != NULL) {
        freeAVLTree(tree->left);
        freeAVLTree(tree->right);
        struct halfEdgeAVL* current = tree;
        while (current != NULL) {
            struct halfEdgeAVL* temp = current;
            current = current->equal;

            free(temp->data->error);
            free(temp->data->vbar);
            free(temp->data);
            free(temp);
        }
    }
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
            current->error = malloc(sizeof(struct halfEdgeError));
            current->vbar = calloc(3, sizeof(float));

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

void computeHalfEdgeMatrixError(struct halfEdge* he) {
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
        he->error->err = error(Qbar, vbar);
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

            he->error->err = err1;
        } else if (err2 < err_mid) {
            he->vbar[0] = v2[0];
            he->vbar[1] = v2[1];
            he->vbar[2] = v2[2];

            he->error->err = err2;
        } else {
            he->vbar[0] = v_mid[0];
            he->vbar[1] = v_mid[1];
            he->vbar[2] = v_mid[2];

            he->error->err = err_mid;
        }
    }
}

float edgeLengthSquared(struct halfEdge* edge) {
    float x_len = edge->endVertex->x - edge->startVertex->x;
    float y_len = edge->endVertex->y - edge->startVertex->y;
    float z_len = edge->endVertex->z - edge->startVertex->z;
    return x_len * x_len + y_len * y_len + z_len * z_len;
}

bool validContract(struct halfEdge* he) {
    struct vertice** oneRingVertices = malloc(sizeof(struct vertice*) * nbVertices);

    oneRingVertices[0] = he->startVertex;

    int index = 1;
    int nbVois = 0;
    struct halfEdge* cur = he->next->pair;
    while (cur != he) {
        oneRingVertices[index] = cur->startVertex;
        cur = cur->next->pair;
        index++;
        nbVois++;
    }

    int cpt = 0;

    for (int i = 0; i < index; i++) {
        if (he->endVertex == oneRingVertices[i]) {
            cpt++;
        }
    }

    cur = he->pair->next;
    while (cur != he) {
        bool seen = false;
        for (int i = 0; i < index; i++) {
            if (cur->endVertex == oneRingVertices[i]) {
                seen = true;
                cpt++;
            }
        }
        if (!seen) {
            nbVois++;
        }
        cur = cur->pair->next;
    }

    free(oneRingVertices);
    return (nbVois != 2) && (cpt == 2);
}

bool goodContract(struct halfEdge* he) {
    struct halfEdge* current = he->next->pair;
    float a, b, c, d;
    float abar, bbar, cbar, dbar;
    struct vertice* A;
    struct vertice* B;
    struct vertice* Bbar = malloc(sizeof(struct vertice));
    struct vertice* C;
    while (current != he->pair->prev) {
        A = he->startVertex;
        B = he->endVertex;
        C = he->next->endVertex;
        computePlaneEquation(A, B, C, &a, &b, &c, &d);

        A = he->startVertex;
        Bbar->x = he->vbar[0];
        Bbar->y = he->vbar[1];
        Bbar->z = he->vbar[2];
        C = he->next->endVertex;
        computePlaneEquation(A, Bbar, C, &abar, &bbar, &cbar, &dbar);

        float dot = a*abar + b*bbar + c*cbar;
        if (dot < 0.8) {
            free(Bbar);
            return false;
        }

        current = current->next->pair;
    }

    current = he->pair->next->pair;
    while (current != he->prev) {
        A = he->startVertex;
        B = he->endVertex;
        C = he->next->endVertex;
        computePlaneEquation(A, B, C, &a, &b, &c, &d);

        A = he->startVertex;
        Bbar->x = he->vbar[0];
        Bbar->y = he->vbar[1];
        Bbar->z = he->vbar[2];
        C = he->next->endVertex;
        computePlaneEquation(A, Bbar, C, &abar, &bbar, &cbar, &dbar);

        float dot = a*abar + b*bbar + c*cbar;
        if (dot < 0.8) {
            free(Bbar);
            return false;
        }

        current = current->next->pair;
    }

    free(Bbar);

    return true;
}

void computeHalfEdgeError(struct halfEdge* he) {
    computeHalfEdgeMatrixError(he);
    he->error->valid = validContract(he);
    he->error->good = goodContract(he);
}

void computeError() {
    struct halfEdgeList* current;
    while (edges != NULL) {

        current = edges;

        computeHalfEdgeError(current->current);

        edgesTree = insertNode(edgesTree, current->current);

        edges = edges->next;
        free(current);
    }
}

void clearGlobals() {
    for (int i = 0; i < verticesLen; i++) {
        if (vertices[i] != NULL) {
            free(vertices[i]->Q);
            free(vertices[i]);
        }
    }
    free(vertices);

    for (int i = 0; i < facesLen; i++) {
        if (faces[i] != NULL) {
            free(faces[i]->vertexIndex);
            free(faces[i]);
        }
    }
    free(faces);

    // struct halfEdgeList* current;
    // while (edges != NULL) {
    //     current = edges;
    //     edges = edges->next;
    //     free(current->current->error);
    //     free(current->current->vbar);
    //     free(current->current);
    //     free(current);
    // }

    freeAVLTree(edgesTree);
}

void computeVertexMatrix(struct halfEdge* he) {
    for (int i = 0; i < 16; i++) he->endVertex->Q[i] = 0;

    float a, b, c, d;
    struct vertice* A = he->startVertex;
    struct vertice* B = he->endVertex;
    struct vertice* C = he->next->endVertex;
    computePlaneEquation(A, B, C, &a, &b, &c, &d);
    addQmatrix(he->endVertex, a, b, c, d);

    struct halfEdge* current = he->next->pair;
    while (current != he) {
        A = current->startVertex;
        B = current->endVertex;
        C = current->next->endVertex;
        computePlaneEquation(A, B, C, &a, &b, &c, &d);
        addQmatrix(he->endVertex, a, b, c, d);
    
        current = current->next->pair;
    }
}

void printHalfEdgeVois(struct halfEdge* he) {
    printf("\n");
    struct halfEdge* cur = he->next->pair->next;
    printf("\t");
    while (cur != he->pair->prev->pair) {
        printf("%3i ", cur->endVertex->index);
        cur = cur->pair->next;
    }
    printf("\n");
    printf("\t%3i     %3i\n", he->next->endVertex->index, he->pair->prev->startVertex->index);
    printf("\t    %3i\n", he->endVertex->index);
    printf("\t     ^ \n");
    printf("\t     | \n");
    printf("\t     | \n");
    printf("\t    %3i\n", he->startVertex->index);
    printf("\t%3i     %3i\n", he->prev->startVertex->index, he->pair->next->endVertex->index);
    cur = he->prev->pair->prev;
    printf("\t");
    while (cur != he->pair->next->pair) {
        printf("%3i ", cur->startVertex->index);
        cur = cur->pair->prev;
    }
    printf("\n");
    printf("\n");
}

void printHalfEdgeStartVertexVois(struct halfEdge* he) {
    struct halfEdge* cur = he->pair->next;
    printf("vertex %3i vois : %3i ", he->startVertex->index, he->endVertex->index);
    while (cur != he) {
        printf("%3i ", cur->endVertex->index);
        cur = cur->pair->next;
    }
    printf("\n\n");
}

void printTreeInOrder(struct halfEdgeAVL* tree) {
    if (tree == NULL) return;
    printTreeInOrder(tree->left);
    struct halfEdgeAVL* curr = tree;
    while (curr != NULL) {
        printf("%s %s %lf\n", curr->data->error->valid ? "t" : "f", curr->data->error->good ? "t" : "f",curr->data->error->err);
        curr = curr->equal;
    }
    printf("\n");
    printTreeInOrder(tree->right);
}

void contractEdge(struct halfEdge* he) {
    struct halfEdge* pair = he->pair;

    edgesTree = deleteNode(edgesTree, findNode(edgesTree, he));
    edgesTree = deleteNode(edgesTree, findNode(edgesTree, pair));

    // ajustement des halfedge au point d'arrive
    struct halfEdge* cur = he->next;
    while (cur != pair) {
        cur->startVertex = he->startVertex;
        cur->pair->endVertex = he->startVertex;
        cur = cur->pair->next;
    }

    // ajustemet des face au point d'arrive
    cur = he->next->pair->next;
    while (cur != pair) {
        for (int j = 0; j < cur->face->len; j++) {
            if (cur->face->vertexIndex[j] == he->endVertex) {
                cur->face->vertexIndex[j] = he->startVertex;
            }  
        }
        cur = cur->pair->next;
    }

    // ajustemet de la face directe
    struct face* heFace = he->face;
    if (heFace->len == 3) {
        he->prev->pair->pair = he->next->pair;
        he->next->pair->pair = he->prev->pair;

        edgesTree = deleteNode(edgesTree, findNode(edgesTree, he->prev));
        edgesTree = deleteNode(edgesTree, findNode(edgesTree, he->next));

        free(he->prev->error);
        free(he->prev->vbar);
        free(he->prev);
        free(he->next->error);
        free(he->next->vbar);
        free(he->next);

        nbFaces -= 1;
        faces[he->face->index] = NULL;
        free(he->face->vertexIndex);
        free(he->face);

    } else {
        he->prev->next = he->next;
        he->next->prev = he->prev;
        
        int j_cop = 0;
        struct vertice** vertexIndexList = calloc(heFace->len - 1, sizeof(struct vertice*));
        for (int j_ori = 0; j_ori < heFace->len; j_ori++) {
            if (heFace->vertexIndex[j_ori] != he->endVertex) {
                vertexIndexList[j_cop] = heFace->vertexIndex[j_ori];
                j_cop++;
            }
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
        edgesTree = deleteNode(edgesTree, findNode(edgesTree, pair->prev));
        edgesTree = deleteNode(edgesTree, findNode(edgesTree, pair->next));

        free(pair->prev->error);
        free(pair->prev->vbar);
        free(pair->prev);
        free(pair->next->error);
        free(pair->next->vbar);
        free(pair->next);
        
        nbFaces -= 1;
        faces[pair->face->index] = NULL;
        free(pair->face->vertexIndex);
        free(pair->face);

    } else {
        pair->prev->next = pair->next;
        pair->next->prev = pair->prev;

        int j_cop = 0;
        struct vertice** vertexIndexList = calloc(pairFace->len - 1, sizeof(struct vertice*));
        for (int j_ori = 0; j_ori < pairFace->len; j_ori++) {
            if (pairFace->vertexIndex[j_ori] != he->endVertex) {
                vertexIndexList[j_cop] = pairFace->vertexIndex[j_ori];
                j_cop++;
            }
        }
        pairFace->len--;
        free(pairFace->vertexIndex);
        pairFace->vertexIndex = vertexIndexList;
    }

    // deplacement du sommet
    he->startVertex->x = he->vbar[0];
    he->startVertex->y = he->vbar[1];
    he->startVertex->z = he->vbar[2];
    computeVertexMatrix(he->prev->pair->pair);


    // suppression du sommet de la liste des sommets
    nbVertices -= 1;
    vertices[he->endVertex->index] = NULL;
    free(he->endVertex->Q);
    free(he->endVertex);

    free(he->error);
    free(pair->error);
    free(he->vbar);
    free(pair->vbar);
    free(pair);
    free(he);
}

struct halfEdgeAVL* findMinErrorEdge(struct halfEdgeAVL* tree) {
    struct halfEdgeAVL* node = minValueNode(tree);
    return node;
}

void recomputeErrorVertex(struct halfEdge* sentinelEdgeForErrorRecompute) {
    struct halfEdge* current = sentinelEdgeForErrorRecompute->next->pair;
    struct halfEdgeAVL* treeNode = findNode(edgesTree, sentinelEdgeForErrorRecompute);
    struct halfEdgeError* oldErr = malloc(sizeof(struct halfEdgeError));
    *oldErr = *(treeNode->data->error);
    computeHalfEdgeError(sentinelEdgeForErrorRecompute);
    updateTreeFromNode(treeNode, oldErr);

    treeNode = findNode(edgesTree, sentinelEdgeForErrorRecompute->pair);
    *oldErr = *(treeNode->data->error);
    sentinelEdgeForErrorRecompute->pair->error->valid = sentinelEdgeForErrorRecompute->error->valid;
    sentinelEdgeForErrorRecompute->pair->error->good = sentinelEdgeForErrorRecompute->error->good;
    sentinelEdgeForErrorRecompute->pair->error->err = sentinelEdgeForErrorRecompute->error->err;
    sentinelEdgeForErrorRecompute->pair->vbar[0] = sentinelEdgeForErrorRecompute->vbar[0];
    sentinelEdgeForErrorRecompute->pair->vbar[1] = sentinelEdgeForErrorRecompute->vbar[1];
    sentinelEdgeForErrorRecompute->pair->vbar[2] = sentinelEdgeForErrorRecompute->vbar[2];
    updateTreeFromNode(treeNode, oldErr);

    while (current != sentinelEdgeForErrorRecompute) {
        treeNode = findNode(edgesTree, current);
        *oldErr = *(treeNode->data->error);
        computeHalfEdgeError(current);
        updateTreeFromNode(treeNode, oldErr);

        treeNode = findNode(edgesTree, current->pair);
        *oldErr = *(treeNode->data->error);
        current->pair->error->valid = current->error->valid;
        current->pair->error->good = current->error->good;
        current->pair->error->err = current->error->err;
        current->pair->vbar[0] = current->vbar[0];
        current->pair->vbar[1] = current->vbar[1];
        current->pair->vbar[2] = current->vbar[2];
        updateTreeFromNode(treeNode, oldErr);

        current = current->next->pair;
    }
    free(oldErr);
}

struct halfEdge* lastMin;

int contractMinErrorEdge() {
    struct halfEdgeAVL* minErrorEdgeNode = findMinErrorEdge(edgesTree);
    struct halfEdge* minErrorEdge = minErrorEdgeNode->data;

    struct halfEdge* sentinelEdgeForErrorRecompute = minErrorEdge->next->pair;
    
    if (lastMin == minErrorEdge && !minErrorEdge->error->valid) {
        return 2;
    }

    if (!validContract(minErrorEdge)) {
        struct halfEdgeAVL* treeNode = findNode(edgesTree, minErrorEdge);
        struct halfEdgeError* oldErr = malloc(sizeof(struct halfEdgeError));
        *oldErr = *(treeNode->data->error);
        computeHalfEdgeError(minErrorEdge);
        updateTreeFromNode(treeNode, oldErr);

        treeNode = findNode(edgesTree, minErrorEdge->pair);
        *oldErr = *(treeNode->data->error);
        minErrorEdge->pair->error->valid = minErrorEdge->error->valid;
        minErrorEdge->pair->error->good = minErrorEdge->error->good;
        minErrorEdge->pair->error->err = minErrorEdge->error->err;
        minErrorEdge->pair->vbar[0] = minErrorEdge->vbar[0];
        minErrorEdge->pair->vbar[1] = minErrorEdge->vbar[1];
        minErrorEdge->pair->vbar[2] = minErrorEdge->vbar[2];
        updateTreeFromNode(treeNode, oldErr);

        free(oldErr);
        lastMin = minErrorEdge;
        return 1;
    }

    if (lastMin != minErrorEdge && !goodContract(minErrorEdge)) {
        struct halfEdgeAVL* treeNode = findNode(edgesTree, minErrorEdge);
        struct halfEdgeError* oldErr = malloc(sizeof(struct halfEdgeError));
        *oldErr = *(treeNode->data->error);
        computeHalfEdgeError(minErrorEdge);
        updateTreeFromNode(treeNode, oldErr);

        treeNode = findNode(edgesTree, minErrorEdge->pair);
        *oldErr = *(treeNode->data->error);
        minErrorEdge->pair->error->valid = minErrorEdge->error->valid;
        minErrorEdge->pair->error->good = minErrorEdge->error->good;
        minErrorEdge->pair->error->err = minErrorEdge->error->err;
        minErrorEdge->pair->vbar[0] = minErrorEdge->vbar[0];
        minErrorEdge->pair->vbar[1] = minErrorEdge->vbar[1];
        minErrorEdge->pair->vbar[2] = minErrorEdge->vbar[2];
        updateTreeFromNode(treeNode, oldErr);

        free(oldErr);
        lastMin = minErrorEdge;
        return 1;
    }

    lastMin = NULL;

    contractEdge(minErrorEdge);
    recomputeErrorVertex(sentinelEdgeForErrorRecompute);
    return 0;
}

void contractEdgeTo(int goalVertices) {
    int nbIterations = nbVertices - goalVertices;
    printf("Contracting edge :\n");
    for (int i = 0; i < nbIterations; i++) {
        printLoadingBar(i, nbIterations);
        
        // 0 contract succes
        // 1 contract to do
        // 2 no contract valid
        bool loop = true;
        while (loop) {
            bool succes = contractMinErrorEdge();
            if (succes == 0) {
                loop = false;
            } else if (succes == 1) {
                loop = true;
            } else if (succes == 2) {
                printf("\n Error : No more edge can be contracted\n");
                return;
            }
        }
        
    }
    printLoadingBar(nbIterations, nbIterations);
    printf("\n\n");
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        // Print usage information if the number of arguments is incorrect
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *fileFrom = argv[1];
    const char *fileDest = argv[2];

    // const char *fileFrom = "3D-Models/eight.off";
    // const char *fileDest = "Simplified/eight_res_ter.off";

    printf("%li, %li, %li\n", sizeof(struct vertice), sizeof(struct face), sizeof(struct halfEdge));
    printf("%li, %li\n", 16*sizeof(float), 3*sizeof(struct vertice*));

    readFileOFF(fileFrom);
    printf("\n");

    int goalVertices;
    printf("Enter the number of vertices for the output: ");
    scanf("%d", &goalVertices);
    printf("\n");
    
    setHalfEdgeStruct();

    computeError();

    // printTreeInOrder(edgesTree);

    contractEdgeTo(goalVertices);

    writeFileOFF(fileDest);

    clearGlobals();
    return 0;
}
