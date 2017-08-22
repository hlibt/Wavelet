#include <iostream>
#include <stdio.h>
using namespace std;

struct node {
    double data;
    node* left;
    node* right;
};

node* get_new_node(double data){
    node* new_node = new node();
    new_node->data = data;
    new_node->left = NULL;
    new_node->right = NULL;
    return new_node;
};

node* insert(node* root, double data) {
    if (root==NULL) {
        root = get_new_node(data);
    } else if (data <= root->data) {
        root->left = insert(root->left,data);
    } else if (data > root->data) {
        root->right = insert(root->right,data);
    }
    return root;
};

void print_tree(struct node* node) {
    if (node==NULL) return;
    print_tree(node->left);
    print_tree(node->right);
    cout << node->data << endl;
 //    printf("%f",node->data);
};

bool lookup(node* root,int data) {
    if (root == NULL) return false;
    else if (root->data == data) return true;
    else if (data <= root->data) return lookup(root->left,data);
    else return lookup(root->right,data);
};

int main() {
    node* root = NULL;          // setting tree as empty
    root = insert(root, -1);
    root = insert(root,-1);      // return new node called root, which contains data point of 4
    root = insert(root,-0.5);
    root = insert(root,-1);
    root = insert(root,-0.75);    
    root = insert(root,-0.5);
    root = insert(root,-0.25);
    print_tree(root);
    return 0;
}
    
