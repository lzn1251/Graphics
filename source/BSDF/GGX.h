#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
class GGX
{
public:
	GGX(const glm::vec2& roughness, const glm::vec3& worldNormal = glm::vec3(0.0f, 0.0f, 1.0f), const glm::vec3& worldTangent = glm::vec3(1.0f, 0.0f, 0.0f))
		: m_Roughness(roughness), m_WorldNormal(worldNormal), m_WorldTangent(worldTangent) {
		m_WorldBitangent = glm::cross(m_WorldNormal, m_WorldTangent);
		m_LocalToWorld = glm::mat3(m_WorldTangent, m_WorldBitangent, m_WorldNormal);
		m_WorldToLocal = glm::inverse(m_LocalToWorld);
	}

	GGX(const GGX&) = delete;
	GGX& operator=(const GGX&) = delete;

	bool SampleAnisoGGX(const glm::vec2& sample, const glm::vec3& viewDir, glm::vec3& outDir, float& eval, float& pdf);

	bool SampleAnisoGGXVisible(const glm::vec2& sample, const glm::vec3& viewDir, glm::vec3& outDir, float& eval, float& pdf);

	bool SampleAnisoGGXVisibleSphericalCaps(const glm::vec2& sample, const glm::vec3& viewDir, glm::vec3& outDir, float& eval, float& pdf);

	bool SampleAnisoGGXVisibleSphericalCapsBounded(const glm::vec2& sample, const glm::vec3& viewDir, glm::vec3& outDir, float& eval, float& pdf);

public:
	glm::vec3 SampleAnisoGGXNormal(const glm::vec2& sample, const glm::vec3& localV);

	glm::vec3 SampleAnisoGGXVisibleNormal(const glm::vec2& sample, const glm::vec3& localV);

	glm::vec3 SampleAnisoGGXVisibleNormalSphericalCaps(const glm::vec2& sample, const glm::vec3& localV);

	glm::vec3 SampleAnisoGGXVisibleNormalSphericalCapsBounded(const glm::vec2& sample, const glm::vec3& localV);

private:
	float DAnisoGGX(const glm::vec3& normal);

	float SmithAnisoGGXShadowG1(const glm::vec3& viewDir);

	float SmithAnisoGGXShadowG2(const glm::vec3& wi, const glm::vec3& wo);

	float Lambda(const glm::vec3& viewDir);

	float GGXBoundedVndfPDF(const glm::vec3& wi, const glm::vec3& wo);

private:
	glm::vec2 m_Roughness;
	glm::vec3 m_WorldNormal;
	glm::vec3 m_WorldTangent;
	glm::vec3 m_WorldBitangent;

	glm::mat3 m_LocalToWorld;
	glm::mat3 m_WorldToLocal;
};
