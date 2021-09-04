package ds.finalProject;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

@SuppressWarnings("unchecked")
public class HashTable<K, V> implements Iterable<HashTable.Node<K, V>> {
    /**
     * For the sake of separate chaining structure, we need a container to hold the (key, value) pairs and
     * a reference to the next connected node (chaining happens here).
     * @param <K> Type of the key
     * @param <V> Type of the value
     */
    static class Node<K, V> {
        private final K key;
        private V value;
        private Node<K, V> next; //Reference to the next connected node (Imagine a chain structure).

        public Node(K key, V value) {
            this.key = key;
            this.value = value;
        }

        public K getKey() {
            return key;
        }

        public V getValue() {
            return value;
        }

        public void setValue(V value) {
            this.value = value;
        }

        public Node<K, V> getNext() {
            return next;
        }

        public void setNext(Node<K, V> next) {
            this.next = next;
        }
    }

    private Node<K, V>[] table; //Key's will be hashed and inserted in this array.
    private int occupied; //Number of occupied spaces of the table.

    public HashTable() {
        table = new Node[16]; //Initial size of the table is set to 16.
        occupied = 0;
    }

    public int getOccupied() {
        return occupied;
    }

    public Node<K, V>[] getTable() {
        return table;
    }

    //Generates a hash code for a string.
    private int hashString(String string) {
        int hash = 7;

        for (int i = 0; i < string.length(); i++) {
            hash = (hash << 5) - hash + string.charAt(i);
        }

        return hash;
    }

    //Uses java's hashCode method to hash the given key. Then, using modulo (%) we obtain the index of the table array
    //in which the given key should be chained.
    private int getTableIndex(K key) {
        int hash = key.hashCode();

        //Because the generated hash might be a negative number and modulo of a negative number will result in a negative
        //index, so, we need to make sure that the hash is a positive number before we get the modulo.
        if (hash < 0) {
            hash *= -1;
        }

        return hash % table.length;
    }

    /**
     * Removes the node containing the given key from the hashTable.
     * @param key the key of the node to be removed.
     * @return the value of the removed node.
     */
    public V remove(K key) {
        int tableIndex = getTableIndex(key);

        Node<K, V> currentNode = table[tableIndex];
        Node<K, V> previousNode = null;

        //Traverses and finds the node containing the given key.
        while (currentNode != null) {
            if (currentNode.getKey().equals(key)) {
                break;
            }

            previousNode = currentNode;
            currentNode = currentNode.getNext();
        }

        //If we couldn't find the node containing the given key.
        if (currentNode == null) {
            throw new IllegalArgumentException("Couldn't find the key: " + key);
        }

        if (previousNode == null) {
            //If the node to be removed is the head of the chain, then, we only need to update the head.
            table[tableIndex] = currentNode.next;
        } else {
            previousNode.next = currentNode.next;
        }

        occupied--;

        return currentNode.getValue();
    }

    /**
     * Returns the value of the node containing the given key.
     * @param key the key of the node to look for.
     * @return the value of the node containing the given key.
     */
    public V get(K key) {
        int tableIndex = getTableIndex(key);

        Node<K, V> currentNode = table[tableIndex];

        //Traverses and finds the node containing the given key. If found, returns the value stored in that node.
        while (currentNode != null) {
            if (currentNode.getKey().equals(key)) {
                return currentNode.getValue();
            }

            currentNode = currentNode.getNext();
        }

        throw new IllegalArgumentException("Couldn't find the key." + key);
    }

    /**
     * Puts a new (key, value) pair in the hashTable. At first it finds the proper index, then searches the chain in that index
     * to see if we already have a node containing the given key. If so, we update the value stored in that node. If a node containing
     * the given key doesn't exist, we construct a new node containing the given key and value pair and make it the head of that chain in the
     * proper index.
     * @param key the key of the new node
     * @param value the value of the node
     */
    public void put(K key, V value) {
        int tableIndex = getTableIndex(key);

        Node<K, V> currentNode = table[tableIndex];

        //Traverses the chain to see if we already have a node containing the given key. If so, update that node's value.
        while (currentNode != null) {
            if (currentNode.getKey().equals(key)) {
                currentNode.setValue(value);
                return;
            }

            currentNode = currentNode.getNext();
        }

        //Constructing a new node containing the given (key, value) pair and making that node the head of the corresponding chain.
        currentNode = table[tableIndex];
        Node<K, V> newNode = new Node<>(key, value);
        newNode.setNext(currentNode);
        table[tableIndex] = newNode;
        occupied++;

        //To keep the chains smaller as possible, we need to handle the load factor.
        //If 70% of the table is filled, we need to extend the size of the table and re hash all the existing keys and store them back in the table.
        if (occupied >= (7 * table.length) / 10) {
            Node<K, V>[] temp = table.clone(); //To keep the previous table for re hashing process.
            table = new Node[2 * table.length];
            occupied = 0;

            //For each chain and for each node of a chain, we put that (key, value) pair again in the table (re hash and put).
            for (Node<K, V> node : temp) {
                while (node != null) {
                    put(node.getKey(), node.getValue());
                    node = node.getNext();
                }
            }
        }
    }

    /**
     * Traverses the table and chains to see if we have a node having the given key.
     * @param key the key to look for.
     * @return true, if the hashTable has a node containing the given key.
     */
    public boolean containsKey(K key) {
        for (Node<K, V> currentNode : table) {
            while (currentNode != null) {
                if (currentNode.getKey().equals(key)) {
                    return true;
                }

                currentNode = currentNode.getNext();
            }
        }

        return false;
    }

    /**
     * Traverses the hashTable and adds all the existing keys to a set.
     * @return a set containing all the existing keys
     */
    public Set<K> keySet() {
        Set<K> set = new HashSet<>();

        for (Node<K, V> currentNode : table) {
            if (currentNode != null) {
                Node<K, V> currentHeadNode = currentNode;
                while (currentHeadNode != null) {
                    set.add(currentNode.getKey());
                    currentHeadNode = currentHeadNode.getNext();
                }
            }
        }

        return set;
    }

    /**
     * Iterates over the proper chain of the table and looks for the node having the given key and return it.
     * @param key the key to look for
     * @return the node containing the given key or it throws an exception if it couldn't find the node.
     */
    public Node<K, V> getNode(K key) {
        int tableIndex = getTableIndex(key);
        Node<K, V> currentNode = table[tableIndex];

        while (currentNode != null) {
            if (currentNode.getKey().equals(key)) {
                return currentNode;
            }

            currentNode = currentNode.getNext();
        }

        throw new IllegalArgumentException("Couldn't find a node having this key: " + key.toString());
    }

    /**
     * @return an iterator of the nodes of the hashTable.
     */
    public Iterator<Node<K, V>> iterator() {
        return new Iterator<>() {
            private HashTable.Node<K, V> currentNode = null;
            private int tableIndex = -1;
            int counter = 0;

            @Override
            public boolean hasNext() {
                return counter < occupied;
            }

            private HashTable.Node<K, V> findNextIndex() {
                tableIndex++;
                for (int i = tableIndex; i < getTable().length; i++) {
                    if (getTable()[i] == null) {
                        tableIndex++;
                    } else {
                        break;
                    }
                }

                return getTable()[tableIndex];
            }

            @Override
            public Node<K, V> next() {
                HashTable.Node<K, V> current = currentNode;

                if (currentNode == null) {
                    current = findNextIndex();
                    currentNode = current;
                    counter++;
                    return current;
                } else if (current.getNext() == null) {
                    current = findNextIndex();
                    currentNode = current;
                    counter++;
                    return current;
                } else {
                    current = current.getNext();
                    currentNode = current;
                    counter++;
                    return current;
                }
            }
        };
    }
}