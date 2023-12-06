---
layout: page
title: Lectures
---

<table>
  <tr>
    <th>Title</th>
  </tr>
{% for lecture in site.pages %}
  {% if lecture.layout == 'page' and lecture.element == 'lecture' %}
    {% unless lecture.title contains "Template" %}
      <tr>
      <td nowrap><a href="{{ lecture.url | prepend: site.baseurl }}">
      {{ lecture.title }}</a></td>
      </tr>
    {% endunless %}
  {% endif %}
{% endfor %}
</table>