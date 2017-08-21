#include <iostream>
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
        return root;
    } else if (data <= root->data) {
        root->left = insert(root->left,data);
    } else if (data > root->data) {
        root->right = insert(root->right,data);
    }
    return root;
};

bool lookup(node* root,int data) {
    if (root == NULL) return false;
    else if (root->data == data) return true;
    else if (data <= root->data) return lookup(root->left,data);
    else return lookup(root->right,data);
};

int main() {
    node* root=NULL;    // setting tree as empty
    insert(root,-1);
    insert(root,-0.75);
       
}
    
