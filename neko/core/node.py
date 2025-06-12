class Node:
    """
    Represents a node in the biological network (gene, protein, metabolite, etc.).
    """
    def __init__(self, node_id: str, node_type: str = "gene", metadata: dict = None):
        self.id = node_id
        self.type = node_type
        self.metadata = metadata or {}

    def __repr__(self):
        return f"Node(id={self.id}, type={self.type}, metadata={self.metadata})"

    def __eq__(self, other):
        return isinstance(other, Node) and self.id == other.id

    def __hash__(self):
        return hash(self.id)

