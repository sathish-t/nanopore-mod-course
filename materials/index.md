---
layout: page
title: Hands-on Materials
---

<table>
  <tr>
    <th>Title</th>
  </tr>
{% for material in site.pages %}
  {% if material.layout == 'page' and material.element == 'notes' %}
    {% unless material.title contains "Template" %}
      <tr>
      <td nowrap><a href="{{ material.url | prepend: site.baseurl }}">
      {{ material.title }}</a></td>
      </tr>
    {% endunless %}
  {% endif %}
{% endfor %}
</table>