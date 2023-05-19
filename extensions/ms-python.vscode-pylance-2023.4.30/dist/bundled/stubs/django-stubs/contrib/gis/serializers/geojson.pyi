from typing import Any

from django.core.serializers.json import Serializer as JSONSerializer

class Serializer(JSONSerializer):
    def start_serialization(self) -> None: ...
    def end_serialization(self) -> None: ...
    geometry_field: Any = ...
    def start_object(self, obj: Any) -> None: ...
    def get_dump_object(self, obj: Any) -> Any: ...
    def handle_field(self, obj: Any, field: Any) -> None: ...

class Deserializer:
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
