class Edge:
    """
    Represents an edge (interaction) in the biological network.
    Stores source, target, interaction type, evidence, database of origin, and metadata.
    """
    def __init__(self, source: str, target: str, interaction_type: str = "undefined",
                 evidence: str = None, database: str = None, metadata: dict = None):
        self.source = source
        self.target = target
        self.interaction_type = interaction_type
        self.evidence = evidence
        self.database = database
        self.metadata = metadata or {}

    def __repr__(self):
        return (f"Edge(source={self.source}, target={self.target}, type={self.interaction_type}, "
                f"evidence={self.evidence}, database={self.database}, metadata={self.metadata})")

    def __eq__(self, other):
        return (isinstance(other, Edge) and self.source == other.source and self.target == other.target
                and self.interaction_type == other.interaction_type)

    def __hash__(self):
        return hash((self.source, self.target, self.interaction_type))

