from __future__ import annotations

from typing import Any

import httpx
from openai import OpenAI


def make_openai_client(api_key: str, base_url: str, timeout: float | None = None) -> OpenAI:
    """Create an OpenAI-compatible client while keeping the project's SSL behavior."""
    http_client = httpx.Client(verify=False, timeout=timeout) if timeout else httpx.Client(verify=False)
    return OpenAI(api_key=api_key, base_url=base_url, http_client=http_client)


def _as_dict(value: Any) -> dict[str, Any] | None:
    if isinstance(value, dict):
        return value
    model_dump = getattr(value, "model_dump", None)
    if callable(model_dump):
        try:
            dumped = model_dump()
            return dumped if isinstance(dumped, dict) else None
        except Exception:
            return None
    return None


def _content_to_text(content: Any) -> str:
    if isinstance(content, str):
        return content.strip()
    if isinstance(content, list):
        parts: list[str] = []
        for item in content:
            if isinstance(item, str):
                parts.append(item)
                continue
            item_dict = _as_dict(item)
            if item_dict:
                text = item_dict.get("text") or item_dict.get("content")
                nested = _content_to_text(text)
                if nested:
                    parts.append(nested)
        return "\n".join(part.strip() for part in parts if part and part.strip()).strip()
    content_dict = _as_dict(content)
    if content_dict:
        for key in ("text", "content", "output_text"):
            nested = _content_to_text(content_dict.get(key))
            if nested:
                return nested
        parts = content_dict.get("parts")
        nested = _content_to_text(parts)
        if nested:
            return nested
    return ""


def extract_chat_text(response: Any) -> str:
    """Extract text from common OpenAI-compatible response variants."""
    try:
        choices = getattr(response, "choices", None) or []
        if choices:
            choice = choices[0]
            message = getattr(choice, "message", None)
            text = _content_to_text(getattr(message, "content", None) if message is not None else None)
            if text:
                return text
            text = _content_to_text(getattr(choice, "text", None))
            if text:
                return text

        response_dict = _as_dict(response) or {}
        text = _content_to_text(response_dict.get("output_text"))
        if text:
            return text

        dict_choices = response_dict.get("choices") or []
        if dict_choices:
            choice_dict = _as_dict(dict_choices[0]) or {}
            message_dict = _as_dict(choice_dict.get("message")) or {}
            text = _content_to_text(message_dict.get("content")) or _content_to_text(choice_dict.get("text"))
            if text:
                return text

        # Gemini native and some OpenAI-compatible gateways expose candidates.
        candidates = response_dict.get("candidates") or []
        if candidates:
            candidate_dict = _as_dict(candidates[0]) or {}
            content_dict = _as_dict(candidate_dict.get("content")) or {}
            text = _content_to_text(content_dict.get("parts"))
            if text:
                return text

        output_items = response_dict.get("output") or []
        if output_items:
            return _content_to_text(output_items)
    except Exception:
        return ""
    return ""


def chat_completion_text(
    client: OpenAI,
    model: str,
    messages: list[dict[str, str]],
    temperature: float = 0.1,
    max_tokens: int | None = None,
) -> str:
    """Run chat.completions and return normalized text across compatible providers."""
    kwargs: dict[str, Any] = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
    }
    if max_tokens is not None:
        kwargs["max_tokens"] = max_tokens
    response = client.chat.completions.create(**kwargs)
    return extract_chat_text(response)
